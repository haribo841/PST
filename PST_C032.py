# Python code to compute capacities for given parameters and produce three separate plots (one per Bj).
# This will run in the notebook environment and display/save the images and a results table.
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Parameters requested by the user
C_primary = 8  # pojemność systemu pierwotnego (AUs)
Ci = np.array([C_primary], dtype=int)
Aalt = 4.0  # ruch oferowany do alternatywy (Erlang)
Bj_list = [0.08, 0.10, 0.12]
amin = 0.8
amax = 1.7
astep = 0.1

out_dir = Path("/mnt/data")
out_dir.mkdir(parents=True, exist_ok=True)

def erlangB(A, n):
    if n <= 0:
        return 1.0 if A > 0 else 0.0
    B = 1.0
    for i in range(1, n + 1):
        B = (A * B) / (i + A * B)
    return B

def riordan_overflow_mean_variance(a, m):
    # a: offered traffic to the group; m: capacity (channels)
    if m <= 0:
        return 0.0, 0.0
    B = erlangB(a, m)
    alpha = a * B
    denom = (m + 1 - a + alpha)
    if denom <= 1e-12:
        var = alpha
    else:
        var = alpha * (1.0 - alpha + a / denom)
    return alpha, var

a_values = np.round(np.arange(amin, amax + 1e-9, astep), 10)

# Storage for full results
cols = ["a", "R_total", "Var_total", "Mean_alt_total", "Var_alt_total"]
for Bj in Bj_list:
    cols += [f"Calt_eq_B{Bj}", f"Calt_hay_B{Bj}"]
results = []

for a in a_values:
    # offered traffic to each primary group
    ai = a * Ci.astype(float)
    Ri_list = []
    Var_list = []
    for idx in range(len(Ci)):
        mean_i, var_i = riordan_overflow_mean_variance(ai[idx], Ci[idx])
        Ri_list.append(mean_i)
        Var_list.append(var_i)
    R_total = float(np.sum(Ri_list))
    Var_total = float(np.sum(Var_list))
    Mean_alt_total = Aalt + R_total
    Var_alt_total = Aalt + Var_total

    row = [a, R_total, Var_total, Mean_alt_total, Var_alt_total]
    for Bj in Bj_list:
        # EqChange
        Calt_eq = 0
        if Mean_alt_total <= 0:
            Calt_eq = 0
        else:
            while True:
                B_alt = erlangB(Mean_alt_total, Calt_eq)
                if B_alt <= Bj or Calt_eq > 20000:
                    break
                Calt_eq += 1
        # Fredericks-Hayward
        if Mean_alt_total <= 0:
            Z = 1.0
        else:
            Z = Var_alt_total / Mean_alt_total
            if Z <= 0:
                Z = 1.0
        A_scaled = Mean_alt_total / Z if Z != 0 else Mean_alt_total
        n_star = 0
        if A_scaled <= 0:
            n_star = 0
        else:
            while True:
                B_scaled = erlangB(A_scaled, n_star)
                if B_scaled <= Bj or n_star > 20000:
                    break
                n_star += 1
        Calt_hay = int(math.ceil(n_star * Z))
        row += [Calt_eq, Calt_hay]
    results.append(row)

# Build DataFrame and save as txt and csv
df = pd.DataFrame(results, columns=cols)
txt_path = out_dir / "wyniki_C8_Aalt4_range_a_0.8_1.7.txt"
csv_path = out_dir / "wyniki_C8_Aalt4_range_a_0.8_1.7.csv"
df.to_csv(csv_path, index=False, sep=';')
with open(txt_path, "w", encoding="utf-8") as f:
    f.write(df.to_string(index=False))

# Create one separate figure per Bj (each figure contains two lines: EqChange and Fredericks-Hayward)
plot_paths = []
for Bj in Bj_list:
    idx_eq = cols.index(f"Calt_eq_B{Bj}")
    idx_hay = cols.index(f"Calt_hay_B{Bj}")
    C_eq = df.iloc[:, idx_eq].values
    C_hay = df.iloc[:, idx_hay].values

    plt.figure(figsize=(8,5))
    plt.plot(a_values, C_eq, marker='o', linestyle='-', label=f'EqChange B={Bj}')
    plt.plot(a_values, C_hay, marker='x', linestyle='--', label=f'Fredericks-Hayward B={Bj}')
    plt.xlabel('a (ruch przypadający na 1 jednostkę pojemności)')
    plt.ylabel('Wymagana pojemność systemu alternatywnego Calt')
    plt.title(f'Porównanie metod (C_primary={C_primary}, Aalt={Aalt}, B={Bj})')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    pfile = out_dir / f'porownanie_C8_Aalt4_B{int(Bj*100)}.png'
    plt.savefig(pfile, dpi=200)
    plot_paths.append(str(pfile))
    plt.show()

# Additionally, create combined plot with all Bj curves (two methods per Bj) on single figure
plt.figure(figsize=(10,6))
for j, Bj in enumerate(Bj_list):
    idx_eq = cols.index(f"Calt_eq_B{Bj}")
    idx_hay = cols.index(f"Calt_hay_B{Bj}")
    C_eq = df.iloc[:, idx_eq].values
    C_hay = df.iloc[:, idx_hay].values
    plt.plot(a_values, C_eq, marker='o', linestyle='-', label=f'EqChange B={Bj}')
    plt.plot(a_values, C_hay, marker='x', linestyle='--', label=f'Fredericks-Hayward B={Bj}')

plt.xlabel('a (ruch przypadający na 1 jednostkę pojemności)')
plt.ylabel('Wymagana pojemność systemu alternatywnego Calt')
plt.title(f'Porównanie metod (C_primary={C_primary}, Aalt={Aalt})')
plt.grid(True)
plt.legend(ncol=2)
plt.tight_layout()
combined_path = out_dir / 'porownanie_C8_Aalt4_all_B.png'
plt.savefig(combined_path, dpi=200)
plt.show()
plot_paths.append(str(combined_path))

# Provide paths and a short summary table
plot_paths, str(txt_path), str(csv_path)plot_paths, str(txt_path), str(csv_path)