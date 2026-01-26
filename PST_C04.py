#!/usr/bin/env python3
"""
Program do wyznaczania prawdopodobieństw blokady w systemie alternatywnym (przelewowym)
przy użyciu metody Kaufman-Roberts jako podstawy obliczeń oraz przyjętej interpretacji
wzoru (4) z załącznika.

Uwagi i założenia (wyraźnie proszę sprawdzić ze składem zadania i literaturą):
- Program liczy rozkład zajętości łącz (P(n), n=0..C) metodą Kaufman-Robertsa
  (klasyczna metoda dla systemów wieloprzepustowych z żądaniami o kilku jednostkach zasobu).
- Zgodnie z obrazem wzoru (4) przyjmuję, że prawdopodobieństwo blokady strumienia i-tego
  obliczamy jako sumę P(n) dla stanów, w których pozostaje mniej niż wymaganych jednostek
  zasobu: E_i = sum_{n = C - c + 1}^{C} P(n), gdzie c wybieramy zgodnie z warunkiem
  przedstawionym w obrazie:
    jeżeli Q > C - ci1[i] -> używamy ci0[i]
    w przeciwnym razie        -> używamy ci1[i]
  (dokładna interpretacja tej reguły powinna być skontrolowana przez użytkownika).
- Aby otrzymać P(n) potrzebujemy stałych wartości rozmiarów żądań i obciążeń. Ponieważ
  w zadaniu są dwie wartości rozmiaru (przedprogowa i poprogowa), obliczam rozkład
  dwukrotnie:
    * system_A: wszystkie żądania mają rozmiary = ci0[i] i obciążenia = ai0[i]
    * system_B: wszystkie żądania mają rozmiary = ci1[i] i obciążenia = ai1[i]
  oraz otrzymuję odpowiadające P_A(n) i P_B(n). Następnie dla każdego strumienia i
  wybieram sumę P_A lub P_B zgodnie z regułą powyżej.
  (To przybliżenie: model z prawdziwym progiem wymaga bardziej złożonej, stanowo-zależnej
   formuły; jednak to rozwiązanie jest praktyczne i często stosowane w ćwiczeniach.)

Wejscie (wewnetrzne prompt/parametry):
  amin, amax, astep - wartości ruchu oferowanego przypadającego na jedn. pojemności
  C - pojemność systemu (liczba jednostek zasobu, integer)
  Q - prog (integer)
  m - liczba klas
  ci0 - lista długości (resource units) w obszarze przedprogowym (m elementow, integer)
  ci1 - lista dlugosci w obszarze poprogowym (m elementow, integer)

Wyjscie:
  plik txt zawierajacy naglowek z parametrami oraz tabelę:
    kolumna 1: a (ruch na jednostke pojemnosci)
    kolejne kolumny: prawdopodobienstwa blokady E_1 ... E_m

Jak uruchomić:
  - Uruchomić jako skrypt python: python3 pojemnosc_systemu.py
  - Program poprosi o wprowadzenie danych na wejściu lub może być uruchamiany z opcjami
    w kodzie (przykład w bloku __main__).

"""

from math import isclose
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')
except AttributeError:
    pass  # older Python: rely on environment or use bytes via .buffer


def kaufman_roberts(rho, size, C):
    """
    Oblicza współczynniki g(n) i prawdopodobieństwa P(n) dla n=0..C metodą Kaufman-Roberts.
    rho: lista obciążeń (Erlangs) dla klas [a1, a2, ...]
    size: lista rozmiarów (liczb całkowitych) zasobu wymaganych przez klasy [s1, s2, ...]
    C: pojemność (int)

    Zwraca: P (lista długości C+1) — P[n] = P(N=n)
    """
    m = len(rho)
    # g[0]=1
    g = [0.0] * (C + 1)
    g[0] = 1.0
    # oblicz g(n)
    for n in range(1, C + 1):
        s = 0.0
        for j in range(m):
            sj = size[j]
            if sj <= n:
                s += rho[j] * g[n - sj]
        g[n] = s / n
    G = sum(g)
    if isclose(G, 0.0):
        raise RuntimeError("Normalizing constant is zero (sprawdz dane wejściowe).")
    P = [val / G for val in g]
    return P


def blocking_from_P(P, C, s):
    """
    Oblicza prawdopodobienstwo blokady dla klasy o rozmiarze s,
    na podstawie P(n): B = sum_{n=C-s+1}^{C} P(n)
    """
    start = max(0, C - s + 1)
    return sum(P[start: C + 1])


def generate_range(a_min, a_max, a_step):
    vals = []
    a = a_min
    # zachowaj tolerancje pływającej arytmetyki
    while a <= a_max + 1e-12:
        vals.append(round(a, 12))
        a += a_step
    return vals


def run_calculation(amin, amax, astep, C, Q, m, ci0, ci1, out_filename):
    # Walidacja
    if m <= 0:
        raise ValueError("m musi byc > 0")
    if len(ci0) != m or len(ci1) != m:
        raise ValueError(f"ci0 and ci1 must have length m; got m={m}, len(ci0)={len(ci0)}, len(ci1)={len(ci1)}")
    if any((c <= 0 or not float(c).is_integer()) for c in ci0 + ci1):
        raise ValueError("ci0 i ci1 musza byc dodatnimi integerami")

    a_values = generate_range(amin, amax, astep)

    # Przygotuj naglowek pliku wynikowego
    # open file with explicit encoding
    with open(out_filename, 'w', encoding='utf-8') as f:
        f.write(f"# Pojemnosc systemu C = {C}\n")
        f.write(f"# Prog Q = {Q}\n")
        f.write("# Żądania w obszarze przedprogowym (ci0): ")
        f.write(' '.join(str(int(x)) for x in ci0) + '\n')
        f.write("# Żądania w obszarze poprogowym (ci1):   ")
        f.write(' '.join(str(int(x)) for x in ci1) + '\n')
        f.write("# Kolumny: a ")
        for i in range(1, m + 1):
            f.write(f" E_{i}")
        f.write("\n")

    # Dla kazdej wartosci a policz ai0 i ai1, P_A, P_B i E_i
    results = []
    for a in a_values:
        # obciazenia ai0 i ai1 dla klas i
        ai0 = []
        ai1 = []
        for i in range(m):
            # ai0 = (a * C) / (ci0[i] * m)
            ai0_i = (a * C) / (ci0[i] * m)
            ai1_i = ai0_i * (ci0[i] / ci1[i])
            ai0.append(ai0_i)
            ai1.append(ai1_i)

        # Oblicz P_A (wszystkie klasy maja rozmiar ci0, obciazenia ai0)
        try:
            P_A = kaufman_roberts(ai0, ci0, C)
        except RuntimeError as e:
            raise RuntimeError(f"Blad w obliczaniu P_A dla a={a}: {e}")
        # Oblicz P_B (wszystkie klasy maja rozmiar ci1, obciazenia ai1)
        try:
            P_B = kaufman_roberts(ai1, ci1, C)
        except RuntimeError as e:
            raise RuntimeError(f"Blad w obliczaniu P_B dla a={a}: {e}")

        # Dla kazdej klasy wybierz odpowiednie P do sumowania zgodnie z rownaniem (4)
        E = []
        for i in range(m):
            thresh = C - ci1[i]
            # warunek z obrazka: "dla Q > C - c_{i1} -> uzyj sumy z ci0"
            if Q > thresh:
                # uzyj P_A i ci0
                Bi = blocking_from_P(P_A, C, ci0[i])
            else:
                Bi = blocking_from_P(P_B, C, ci1[i])
            E.append(Bi)

        results.append((a, E))

        # Dopisz do pliku
        with open(out_filename, 'a', encoding='utf-8') as f:
            row = f"{a:.6f}"
            for val in E:
                row += f" {val:.8e}"
            f.write(row + "\n")

    return results


if __name__ == '__main__':
    # Prosty interaktywny tryb użytkownika — mozna zastapic danymi testowymi
    print("Program do wyznaczania prawdopodobienstw blokady (wersja CLI).\n")
    try:
        amin = float(input("amin (minimalny ruch na 1 jednostke pojemnosci): "))
        amax = float(input("amax (maksymalny ruch na 1 jednostke pojemnosci): "))
        astep = float(input("astep (krok obliczen): "))
        C = int(input("C (pojemnosc systemu, integer): "))
        Q = int(input("Q (prog, integer): "))
        print("Podaj ci0 (m liczb calkowitych, rozdzielonych spacja):")
        ci0 = list(map(int, input().strip().split()))
        print("Podaj ci1 (m liczb calkowitych, rozdzielonych spacja):")
        ci1 = list(map(int, input().strip().split()))
        if len(ci0) != len(ci1):
            raise SystemExit("ci0 and ci1 must have same length")
        m = len(ci0)
        out_filename = input("Nazwa pliku wynikowego (np. wynik.txt): ").strip() or 'wynik.txt'

        res = run_calculation(amin, amax, astep, C, Q, m, ci0, ci1, out_filename)
        print(f"Obliczenia zakonczone. Wyniki zapisane w pliku: {out_filename}")
    except Exception as e:
        print("Wystapil blad:", e, file=sys.stderr)
        sys.exit(1)
