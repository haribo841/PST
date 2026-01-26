#ifndef BLOCKING_H
#define BLOCKING_H

#include <vector>
#include <string>
#include <utility>

namespace loss {

    /// Typ wyniku: para (wartość a, wektor E_i dla i=1..m)
    using ResultRow = std::pair<double, std::vector<double>>;

    /// Generuje zakres wartości a od amin do amax z krokiem astep.
    /// Zwraca wektor wartości typu double.
    /// Rzuca std::invalid_argument gdy astep <= 0.
    std::vector<double> generate_range(double amin, double amax, double astep);

    /// Implementacja Kaufman–Roberts:
    /// - A: wektor natężeń (obciążeń) dla klas (rozmiar m)
    /// - r: wektor rozmiarów (wymagań zasobu) dla klas (integery, rozmiar m)
    /// - C: pojemność systemu (integer)
    /// Zwraca rozkład stanu P_n dla n=0..C jako vector<double> rozmiaru C+1.
    /// Może rzucać std::invalid_argument lub std::runtime_error (np. problemy z normalizacją).
    std::vector<double> kaufman_roberts(const std::vector<double>& A,
        const std::vector<int>& r,
        int C);

    /// Oblicza prawdopodobieństwo blokady klasy o rozmiarze r przy znanym rozkładzie P (rozmiar C+1).
    /// Interpretacja: blokada = suma P[n] dla n = C - r + 1 .. C.
    /// Zwraca double w [0,1]. Rzuca std::invalid_argument jeśli rozmiary wektorów są niezgodne.
    double blocking_from_P(const std::vector<double>& P, int C, int r);

    /// Główna funkcja obliczeniowa
    /// Parametry:
    ///   amin, amax, astep - zakres wartości a (double)
    ///   C, Q, m           - parametry całkowite
    ///   ci0, ci1          - wektory int długości m (przed- i poprogowe rozmiary klas)
    ///   out_filename      - nazwa pliku, do którego dopisywane będą wyniki
    /// Zwraca vector<ResultRow> — listę par (a, E_vector).
    /// Rzuca std::invalid_argument lub std::runtime_error w przypadku błędów wejścia/wyjścia lub obliczeń.
    std::vector<ResultRow> run_calculation(double amin, double amax, double astep,
        int C, int Q, int m,
        const std::vector<int>& ci0,
        const std::vector<int>& ci1,
        const std::string& out_filename);

} // namespace loss

#endif // BLOCKING_H
