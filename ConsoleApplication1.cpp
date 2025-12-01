#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "PST_C02.h" // Dołączenie naszego pliku nagłówkowego

using namespace std;

int main() {
    double a_min, a_max, a_step;
    int C, m;

    cout << "Podaj pojemnosc systemu (C): ";
    cin >> C;

    cout << "Podaj liczbe klas strumieni (m): ";
    cin >> m;

    vector<StreamClass> streams(m);
    for (int i = 0; i < m; ++i) {
        cout << "Podaj zadania t_" << (i + 1) << " (liczba AU): ";
        cin >> streams[i].t;
    }

    cout << "Podaj a_min (min ruch na jednostke): ";
    cin >> a_min;

    cout << "Podaj a_max (max ruch na jednostke): ";
    cin >> a_max;

    cout << "Podaj a_step (krok obliczen): ";
    cin >> a_step;

    // Otwarcie plików wynikowych
    ofstream file_block("wyniki_blokada.txt");
    ofstream file_resources("wyniki_zasoby.txt");

    if (!file_block.is_open() || !file_resources.is_open()) {
        cerr << "Blad otwarcia plikow do zapisu!" << endl;
        return 1;
    }

    file_block << fixed << setprecision(6);
    file_resources << fixed << setprecision(5);

    // Nagłówki plików
    file_block << "# C = " << C << "\n#\n";
    for (int i = 0; i < m; ++i) file_block << "# t[" << i + 1 << "] = " << streams[i].t << "\n";
    file_block << "#\na";
    for (int i = 0; i < m; ++i) file_block << "\tt" << i + 1;
    file_block << "\n";

    file_resources << "# C = " << C << "\n#\n";
    for (int i = 0; i < m; ++i) file_resources << "# t[" << i + 1 << "] = " << streams[i].t << "\n";
    file_resources << "#";

    // Główna pętla po wartościach ruchu 'a'
    double current_a = a_min;
    while (current_a <= a_max + 1e-9) {

        // Wyznaczenie a_i
        for (int i = 0; i < m; ++i) {
            streams[i].a_i = (current_a * C) / (m * streams[i].t);
        }

        // Użycie funkcji z PST_C02.cpp
        vector<double> s = calculate_S_distribution(C, m, streams);
        vector<double> p = normalize_probabilities(s);

        // Zapis do pliku blokady
        file_block << current_a;
        for (int i = 0; i < m; ++i) {
            double E_i = 0.0;
            // Suma p[n] dla stanów blokujących (C - t_i + 1 do C)
                int start_idx = C - streams[i].t + 1;
            if (start_idx < 0) start_idx = 0;

            for (int n = start_idx; n <= C; ++n) {
                E_i += p[n];
            }
            file_block << "\t" << E_i;
        }
        file_block << "\n";

        // Zapis do pliku zasobów (średnia liczba zgłoszeń)
        file_resources << "\n" << setprecision(1) << current_a << setprecision(5) << "\n";
        file_resources << "n";
        for (int i = 0; i < m; ++i) file_resources << "\tt" << i + 1;
        file_resources << "\tsum\n";

        for (int n = 0; n <= C; ++n) {
            file_resources << n;
            double row_sum = 0.0;
            for (int i = 0; i < m; ++i) {
                double y_i = 0.0;
                // Obliczanie y_i(n)
                    if ((n - streams[i].t >= 0) && (p[n] > 1e-20)) {
                        y_i = (streams[i].a_i * streams[i].t * p[n - streams[i].t]) / p[n];
                    }
                file_resources << "\t" << y_i;
                row_sum += y_i;
            }
            file_resources << "\t" << row_sum << "\n";
        }

        current_a += a_step;
    }

    file_block.close();
    file_resources.close();

    cout << "Zakonczono. Wyniki w plikach txt." << endl;
    return 0;
}