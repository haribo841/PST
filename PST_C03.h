#ifndef WYMIAROWANIE_HPP
#define WYMIAROWANIE_HPP

#include <vector>
#include <string>
#include <utility>

namespace wymiarowanie {

    /**
     * Domyślne ograniczenie przy wyszukiwaniu liczby kanałów.
     */
    inline constexpr int DEFAULT_MAX_CHANNELS = 20000;

    /**
     * Parametry wejściowe dla obliczeń.
     * - amin, amax, astep : zakres i krok parametrów a
     * - Aalt : ruch kierowany bezpośrednio do alternatywy
     * - Ci : wektor pojemności zasobów pierwotnych (int)
     * - Bj_list : lista dopuszczalnych prawdopodobieństw strat
     * - out_txt, out_png : nazwy plików wyjściowych
     * - max_channels : limit bezpieczeństwa przy poszukiwaniu kanałów
     */
    struct Params {
        double amin = 0.2;
        double amax = 1.3;
        double astep = 0.1;

        double Aalt = 0.0;

        std::vector<int> Ci = { 10, 8, 12 };
        std::vector<double> Bj_list = { 0.01, 0.05 };

        std::string out_txt = "wyniki_pojemnosci_alternatywnych.txt";
        std::string out_png = "porownanie_metod.png";

        int max_channels = DEFAULT_MAX_CHANNELS;
    };

    /**
     * Stabilna rekurencyjna postać Erlang-B:
     * Blocking probability dla obciążenia A i n kanałów.
     * Działa również dla n == 0.
     */
    double erlangB(double A, int n);

    /**
     * Riordan: średnia i wariancja ruchu przelewowego z grupy pierwotnej.
     * Zwraca parę (mean_overflow, var_overflow).
     * Parametry:
     *  - a : offered traffic do grupy (Erlang)
     *  - m : pojemność grupy (int)
     */
    std::pair<double, double> riordan_overflow_mean_variance(double a, int m);

    /**
     * Znajduje najmniejsze n takie, że ErlangB(A, n) <= B_target.
     * Jeśli A == 0 => 0. Ograniczone przez max_channels.
     */
    int minimal_capacity_for_B(double A, double B_target, int max_channels = DEFAULT_MAX_CHANNELS);

    /**
     * Pomocnicze zaokrąglenie wartości typu double (uniknięcie akumulacji błędów
     * przy generowaniu wektora wartości a).
     */
    double round10(double x);

    /**
     * Główna funkcja wykonująca obliczenia:
     *  - iteruje przez wartości 'a' (amin..amax krok astep)
     *  - oblicza Riordan overflow dla każdej jednostki Ci
     *  - porównuje metody (ekwiwalentne zmiany vs Fredericks-Hayward)
     *  - zapisuje wyniki do pliku tekstowego oraz generuje plik danych dla wykresu
     *
     * Implementacja tej funkcji odczytuje parametry z przekazanego obiektu Params.
     */
    void compute_and_write(const Params& params);

} // namespace wymiarowanie

#endif
