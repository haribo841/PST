#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wymiarowanie systemu alternatywnego:
 - metoda ekwiwalentnych zmian (proste sumowanie przelewów -> Erlang-B)
 - metoda Fredericks-Hayward (skalowanie przez peakedness -> Erlang-B)

Wejście: parametry w sekcji 'PARAMETRY'.
Wyjście: plik tekstowy z tabelą wyników oraz plik PNG z wykresem.
"""
from typing import Tuple, List
import math
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# PARAMETRY (dostosuj)
# -------------------------
amin: float = 0.2
amax: float = 1.3
astep: float = 0.1

# ruch kierowany bezpośrednio do zasobów alternatywnych (Erlang)
Aalt: float = 0.0

# liczba zasobów pierwotnych oraz ich pojemności Ci (wektor długości k)
Ci: np.ndarray = np.array([10, 8, 12], dtype=int)
k: int = len(Ci)

# lista dopuszczalnych prawdopodobieństw strat w systemie przelewowym (Bj)
Bj_list: List[float] = [0.01, 0.05]

# pliki wyjściowe
out_txt: str = "wyniki_pojemnosci_alternatywnych.txt"
out_png: str = "porownanie_metod.png"

# zabezpieczenie maksymalnej liczby kanałów przy wyszukiwaniu
MAX_CHANNELS: int = 20000


# -------------------------
# FUNKCJE UŻYTECZNE
# -------------------------
def erlangB(A: float, n: int) -> float:
    """
    Stabilna rekurencyjna postać Erlang-B:
    blocking probability dla obciążenia A i n kanałów.
    Działa również dla n == 0.
    """
    if n <= 0:
        return 1.0 if A > 0.0 else 0.0
    B = 1.0
    # rekurencja B_i = (A * B_{i-1}) / (i + A * B_{i-1}), dla i=1..n
    for i in range(1, n + 1):
        denom = i + A * B
        if denom == 0:
            B = 1.0
        else:
            B = (A * B) / denom
    return float(B)


def riordan_overflow_mean_variance(a: float, m: int) -> Tuple[float, float]:
    """
    Riordan: średnia i wariancja ruchu przelewowego z grupy pierwotnej.
    Parametry: a - offered traffic do grupy (Erlang), m - pojemność (int).
    Zwraca (mean_overflow, var_overflow).
    Uwaga: formuła przyjmuje rozsądne parametry (zwykle a < m+1).
    """
    if m <= 0:
        return 0.0, 0.0
    B = erlangB(a, m)
    alpha = a * B
    denom = (m + 1 - a + alpha)
    if denom <= 1e-12:
        # zabezpieczenie przed dzieleniem przez ~0
        var = alpha
    else:
        var = alpha * (1.0 - alpha + a / denom)
    return float(alpha), float(var)


def minimal_capacity_for_B(A: float, B_target: float, max_channels: int = MAX_CHANNELS) -> int:
    """
    Znajduje najmniejsze n takie, że ErlangB(A, n) <= B_target.
    Jeśli A == 0 => 0. Ograniczone przez max_channels.
    """
    if A <= 0.0:
        return 0
    n = 0
    while True:
        B_val = erlangB(A, n)
        if B_val <= B_target or n >= max_channels:
            break
        n += 1
    return n


# -------------------------
# GŁÓWNA PĘTLA OBLICZEŃ
# -------------------------
def main():
    a_values = np.round(np.arange(amin, amax + 1e-12, astep), 10)
    results: List[List[float]] = []

    for a in a_values:
        # ruch przypadający na każdy zasób pierwotny
        ai = a * Ci.astype(float)

        Ri_list: List[float] = []
        Var_list: List[float] = []
        for idx in range(k):
            mean_i, var_i = riordan_overflow_mean_variance(ai[idx], int(Ci[idx]))
            Ri_list.append(mean_i)
            Var_list.append(var_i)

        R_total = float(np.sum(Ri_list))
        Var_total = float(np.sum(Var_list))
        # dodajemy Aalt (zakładamy Poissona: mean=Aalt, var=Aalt)
        Mean_alt_total = Aalt + R_total
        Var_alt_total = Aalt + Var_total

        row: List[float] = [a]
        for Bj in Bj_list:
            # --- metoda ekwiwalentnych zmian ---
            A_equiv = Mean_alt_total
            Calt_eq = minimal_capacity_for_B(A_equiv, Bj)

            # --- metoda Fredericks-Hayward ---
            if Mean_alt_total <= 0.0:
                Z = 1.0
            else:
                Z = Var_alt_total / Mean_alt_total
                if Z <= 0.0:
                    Z = 1.0

            A_scaled = Mean_alt_total / Z
            n_star = minimal_capacity_for_B(A_scaled, Bj)
            # wynikowa pojemność: zaokrąglenie w górę n_star * Z
            Calt_hay = int(math.ceil(n_star * Z))

            row.append(Calt_eq)
            row.append(Calt_hay)

        results.append(row)

    # -------------------------
    # ZAPIS DO PLIKU TXT
    # -------------------------
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("WYMIAROWANIE SYSTEMU ALTERNATYWNEGO\n")
        f.write("Pojemności zasobów pierwotnych (Ci): " + ", ".join(map(str, Ci.tolist())) + "\n")
        f.write("Ruch kierowany bezpośrednio do alternatywy Aalt = {:.6g}\n".format(Aalt))
        f.write("Dopuszczalne prawdopodobieństwa strat (Bj): " + ", ".join(map(str, Bj_list)) + "\n\n")
        f.write("Tabela wyników:\n")
        hdr = ["a"]
        for Bj in Bj_list:
            hdr.append(f"Calt_eq_B={Bj}")
            hdr.append(f"Calt_hay_B={Bj}")
        f.write("\t".join(hdr) + "\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"Wyniki zapisano do pliku: {out_txt}")

    # -------------------------
    # RYSOWANIE WYKRESU PORÓWNAWCZEGO
    # -------------------------
    plt.figure(figsize=(9, 6))
    for j, Bj in enumerate(Bj_list):
        C_eq = [r[1 + 2 * j] for r in results]
        C_hay = [r[1 + 2 * j + 1] for r in results]
        plt.plot(a_values, C_eq, marker='o', linestyle='-', label=f'EqChange B={Bj}')
        plt.plot(a_values, C_hay, marker='x', linestyle='--', label=f'Fredericks-Hayward B={Bj}')

    plt.xlabel('a (ruch przypadający na 1 jednostkę pojemności)')
    plt.ylabel('Wymagana pojemność systemu alternatywnego Calt')
    plt.title('Porównanie metod: ekwiwalentnych zmian vs Fredericks–Hayward')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"Wykres zapisano do pliku: {out_png}")
    plt.show()


if __name__ == "__main__":
    main()
