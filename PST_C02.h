#ifndef PST_C02_H
#define PST_C02_H

#include <vector>

// Struktura przechowuj¹ca parametry pojedynczego strumienia
struct StreamClass {
    int t;      // Liczba ¿¹danych AU
    double a_i; // Œrednie natê¿enie ruchu oferowanego dla klasy 'i'
};

// Deklaracja funkcji obliczaj¹cej rozk³ad s[n] (Kaufman-Roberts)
std::vector<double> calculate_S_distribution(int C, int m, const std::vector<StreamClass>& streams);

// Deklaracja funkcji normalizuj¹cej s[n] do p[n]
std::vector<double> normalize_probabilities(const std::vector<double>& s);

#endif // PST_C02_H