#include "PST_C02.h"
#include <numeric> // dla std::accumulate

using namespace std;

// Implementacja algorytmu Kaufmana-Robertsa
vector<double> calculate_S_distribution(int C, int m, const vector<StreamClass>& streams) {
    vector<double> s(C + 1, 0.0);
    s[0] = 1.0; // s[0] = 1

    for (int n = 1; n <= C; ++n) {
        double sum = 0.0;
        for (int i = 0; i < m; ++i) {
            // Sprawdzenie warunku n >= t_i, aby nie wyjść poza zakres (p[n-t_i] = 0 dla n < t_i)
            if (n >= streams[i].t) {
                // Wzór (4): s[n] = 1/n * suma(a_i * t_i * s[n-ti])
                sum += streams[i].a_i * streams[i].t * s[n - streams[i].t];
            }
        }
        s[n] = sum / n;
    }
    return s;
}

// Normalizacja wyników, aby suma prawdopodobieństw wynosiła 1
vector<double> normalize_probabilities(const vector<double>& s) {
    double sum_s = accumulate(s.begin(), s.end(), 0.0);
    vector<double> p(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        p[i] = s[i] / sum_s;
    }
    return p;
}