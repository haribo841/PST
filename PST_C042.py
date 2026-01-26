import matplotlib.pyplot as plt

def kaufman_roberts(rho, size, C):
    m = len(rho)
    g = [0.0] * (C + 1)
    g[0] = 1.0
    for n in range(1, C + 1):
        s = 0.0
        for j in range(m):
            if size[j] <= n:
                s += rho[j] * g[n - size[j]]
        g[n] = s / n
    G = sum(g)
    return [x / G for x in g]

def blocking(P, C, c):
    return sum(P[C - c + 1 : C + 1])

# Dane systemu
C = 20
a = 0.8
m = 2
ci0 = [3, 2]
ci1 = [2, 1]

# Natężenia ruchu
ai0 = [(a * C) / (ci0[i] * m) for i in range(m)]
ai1 = [ai0[i] * (ci0[i] / ci1[i]) for i in range(m)]

# Rozkłady stanów
P0 = kaufman_roberts(ai0, ci0, C)
P1 = kaufman_roberts(ai1, ci1, C)

Q_vals = list(range(8, 19))
B1, B2 = [], []

for Q in Q_vals:
    # zgodnie ze wzorem (4)
    if Q > C - ci1[0]:
        B1.append(blocking(P0, C, ci0[0]))
    else:
        B1.append(blocking(P1, C, ci1[0]))

    if Q > C - ci1[1]:
        B2.append(blocking(P0, C, ci0[1]))
    else:
        B2.append(blocking(P1, C, ci1[1]))

# Wykres
plt.figure()
plt.plot(Q_vals, B1, marker='o', label='Klasa 1')
plt.plot(Q_vals, B2, marker='s', label='Klasa 2')
plt.xlabel('Próg Q [AUs]')
plt.ylabel('Prawdopodobieństwo blokady')
plt.title('Wpływ progu Q na blokadę (a = 0.8, C = 20)')
plt.grid(True)
plt.legend()
plt.show()
