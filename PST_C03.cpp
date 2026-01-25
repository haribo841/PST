#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include <cstdlib> // system()
#include <limits>

#include "wymiarowanie.hpp" // deklaracje w namespace wymiarowanie (Params, prototypy)
#include "PST_C03.h"        // (opcjonalnie) do innych typów projektu

namespace wymiarowanie {

    // domyślne ograniczenie
    constexpr int DEFAULT_MAX_CHANNELS = 20000;

    // -------------------------
    // Funkcje numeryczne
    // -------------------------
    double erlangB(double A, int n) {
        if (n <= 0) return (A > 0.0) ? 1.0 : 0.0;
        double B = 1.0;
        for (int i = 1; i <= n; ++i) {
            double denom = i + A * B;
            if (denom == 0.0) {
                B = 1.0;
            }
            else {
                B = (A * B) / denom;
            }
        }
        return B;
    }

    std::pair<double, double> riordan_overflow_mean_variance(double a, int m) {
        if (m <= 0) return { 0.0, 0.0 };
        double B = erlangB(a, m);
        double alpha = a * B;
        double denom = (m + 1.0 - a + alpha);
        double var;
        if (denom <= 1e-12) {
            var = alpha;
        }
        else {
            var = alpha * (1.0 - alpha + a / denom);
        }
        return { alpha, var };
    }

    int minimal_capacity_for_B(double A, double B_target, int max_channels) {
        if (A <= 0.0) return 0;
        int n = 0;
        while (true) {
            double Bval = erlangB(A, n);
            if (Bval <= B_target || n >= max_channels) break;
            ++n;
        }
        return n;
    }

    double round10(double x) {
        return std::round(x * 1e10) / 1e10;
    }

    // -------------------------
    // Główna funkcja: compute_and_write
    // -------------------------
    void compute_and_write(const Params& params) {
        // Pobieramy parametry z obiektu params
        const double amin = params.amin;
        const double amax = params.amax;
        const double astep = params.astep;
        const double Aalt = params.Aalt;
        const std::vector<int>& Ci = params.Ci;
        const std::vector<double>& Bj_list = params.Bj_list;
        const std::string out_txt = params.out_txt;
        const std::string out_png = params.out_png;
        const std::string tmp_dat = params.tmp_dat;
        const std::string tmp_plt = params.tmp_plt;
        const int max_channels = params.max_channels > 0 ? params.max_channels : DEFAULT_MAX_CHANNELS;

        // przygotowanie wektora wartości a
        std::vector<double> a_values;
        for (double a = amin; a <= amax + 1e-12; a += astep) {
            a_values.push_back(round10(a));
        }

        // results: dla każdego a przechowujemy: a, (Calt_eq dla Bj1), (Calt_hay dla Bj1), ...
        std::vector<std::vector<long long>> results;
        results.reserve(a_values.size());

        for (double a : a_values) {
            // 1) ruch przypadający na każdy zasób pierwotny
            std::vector<double> ai;
            ai.reserve(Ci.size());
            for (size_t i = 0; i < Ci.size(); ++i) ai.push_back(a * static_cast<double>(Ci[i]));

            // 2) dla każdego zasobu oblicz overflow mean i variance (Riordan)
            std::vector<double> Ri_list;
            std::vector<double> Var_list;
            Ri_list.reserve(Ci.size());
            Var_list.reserve(Ci.size());

            for (size_t i = 0; i < Ci.size(); ++i) {
                auto pr = riordan_overflow_mean_variance(ai[i], Ci[i]);
                Ri_list.push_back(pr.first);
                Var_list.push_back(pr.second);
            }

            double R_total = 0.0;
            double Var_total = 0.0;
            for (double v : Ri_list) R_total += v;
            for (double v : Var_list) Var_total += v;

            double Mean_alt_total = Aalt + R_total;
            double Var_alt_total = Aalt + Var_total;

            std::vector<long long> row;
            row.push_back(static_cast<long long>(std::llround(a * 10000000000.0))); // opcjonalne scaled a
            for (double Bj : Bj_list) {
                // metoda ekwiwalentnych zmian
                int Calt_eq = minimal_capacity_for_B(Mean_alt_total, Bj, max_channels);

                // metoda Fredericks-Hayward
                double Z;
                if (Mean_alt_total <= 0.0) {
                    Z = 1.0;
                }
                else {
                    Z = Var_alt_total / Mean_alt_total;
                    if (Z <= 0.0) Z = 1.0;
                }
                double A_scaled = (Z == 0.0) ? 0.0 : (Mean_alt_total / Z);
                int n_star = minimal_capacity_for_B(A_scaled, Bj, max_channels);
                long long Calt_hay = static_cast<long long>(std::llround(std::ceil(n_star * Z)));

                row.push_back(static_cast<long long>(Calt_eq));
                row.push_back(Calt_hay);
            }
            results.push_back(row);
        }

        // -------------------------
        // Zapis do pliku tekstowego
        // -------------------------
        std::ofstream fout(out_txt, std::ios::out);
        if (!fout) {
            std::cerr << "Błąd otwarcia pliku do zapisu: " << out_txt << "\n";
            return;
        }

        fout << "WYMIAROWANIE SYSTEMU ALTERNATYWNEGO\n";
        fout << "Pojemnosci zasobow pierwotnych (Ci): ";
        for (size_t i = 0; i < Ci.size(); ++i) {
            fout << Ci[i];
            if (i + 1 < Ci.size()) fout << ", ";
        }
        fout << "\n";
        fout << "Ruch kierowany bezposrednio do alternatywy Aalt = " << Aalt << "\n";
        fout << "Dopuszczalne prawdopodobienstwa strat (Bj): ";
        for (size_t i = 0; i < Bj_list.size(); ++i) {
            fout << Bj_list[i];
            if (i + 1 < Bj_list.size()) fout << ", ";
        }
        fout << "\n\nTabela wynikow:\n";

        // naglowek kolumn
        fout << std::setw(12) << "a";
        for (size_t j = 0; j < Bj_list.size(); ++j) {
            std::ostringstream h1, h2;
            h1 << "Calt_eq_B=" << Bj_list[j];
            h2 << "Calt_hay_B=" << Bj_list[j];
            fout << "\t" << std::setw(12) << h1.str() << "\t" << std::setw(12) << h2.str();
        }
        fout << "\n";

        // wiersze
        for (size_t idx = 0; idx < results.size(); ++idx) {
            double aval = a_values[idx];
            fout << std::fixed << std::setprecision(10) << std::setw(12) << aval;
            for (size_t col = 1; col < results[idx].size(); ++col) {
                fout << "\t" << std::setw(12) << results[idx][col];
            }
            fout << "\n";
        }
        fout.close();
        std::cout << "Wyniki zapisano do pliku: " << out_txt << "\n";

        // -------------------------
        // Przygotowanie pliku danych do wykresu (plot_data.dat)
        // -------------------------
        std::ofstream dfile(tmp_dat, std::ios::out);
        if (!dfile) {
            std::cerr << "Błąd otwarcia pliku do zapisu: " << tmp_dat << "\n";
            return;
        }
        for (size_t i = 0; i < a_values.size(); ++i) {
            dfile << std::fixed << std::setprecision(10) << a_values[i];
            for (size_t col = 1; col < results[i].size(); ++col) {
                dfile << " " << results[i][col];
            }
            dfile << "\n";
        }
        dfile.close();

        // -------------------------
        // Skrypt gnuplot do wygenerowania PNG
        // -------------------------
        std::ofstream plt(tmp_plt, std::ios::out);
        if (!plt) {
            std::cerr << "Błąd otwarcia pliku do zapisu: " << tmp_plt << "\n";
            return;
        }
        plt << "set terminal pngcairo size 900,600 enhanced font 'Arial,10'\n";
        plt << "set output '" << out_png << "'\n";
        plt << "set title 'Porownanie metod: ekwiwalentnych zmian vs Fredericks–Hayward'\n";
        plt << "set xlabel 'a (ruch przypadajacy na 1 jednostke pojemnosci)'\n";
        plt << "set ylabel 'Wymagana pojemnosc systemu alternatywnego Calt'\n";
        plt << "set grid\n";
        plt << "plot ";
        for (size_t j = 0; j < Bj_list.size(); ++j) {
            int col_eq = 2 + 2 * j;
            int col_hay = 2 + 2 * j + 1;
            plt << "'" << tmp_dat << "' using 1:" << col_eq
                << " with linespoints lt 1 pt 7 lw 2 title 'EqChange B=" << Bj_list[j] << "'";
            plt << ", ";
            plt << "'" << tmp_dat << "' using 1:" << col_hay
                << " with linespoints lt 2 pt 4 dashtype 2 lw 2 title 'Fredericks-Hayward B=" << Bj_list[j] << "'";
            if (j + 1 < Bj_list.size()) plt << ", ";
        }
        plt << "\n";
        plt.close();

        int sysres = std::system(("gnuplot " + tmp_plt).c_str());
        if (sysres == 0) {
            std::cout << "Wykres zapisano do pliku: " << out_png << "\n";
        }
        else {
            std::cout << "Nie udało się uruchomić gnuplot (kod: " << sysres << ").\n";
            std::cout << "Plik z danymi: " << tmp_dat << " oraz skrypt: " << tmp_plt << " utworzono.\n";
            std::cout << "Jeśli masz zainstalowany gnuplot, uruchom: gnuplot " << tmp_plt << "\n";
        }
    }

} // namespace wymiarowanie
