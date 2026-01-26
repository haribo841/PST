import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

data_text = """\
0.200000 8.68244848e-46 4.68878855e-46 6.56671992e-47 1.05980474e-47
0.300000 2.03617253e-37 1.17623467e-37 1.98057338e-38 3.52800587e-39
0.400000 2.40545515e-31 1.45943746e-31 2.80770541e-32 5.37298660e-33
0.500000 1.34426312e-26 8.47985949e-27 1.81262849e-27 3.67215692e-28
0.600000 9.92724153e-23 6.46918412e-23 1.50961827e-23 3.20790672e-24
0.700000 1.70048769e-19 1.13963078e-19 2.86812143e-20 6.35254467e-21
0.800000 9.53152724e-17 6.54766433e-17 1.76144550e-17 4.04778728e-18
0.900000 2.19388953e-14 1.54086636e-14 4.40090824e-15 1.04561931e-15
1.000000 2.42961807e-12 1.74113602e-12 5.25134543e-13 1.28646313e-13
1.100000 1.45288913e-10 1.06059651e-10 3.36322767e-11 8.47665613e-12
1.200000 5.11716636e-09 3.79985480e-09 1.26233259e-09 3.26739416e-10
1.300000 1.13536367e-07 8.56600482e-08 2.97210370e-08 7.88858390e-09
"""

rows = []
for line in data_text.strip().splitlines():
    parts = line.split()
    a = float(parts[0])
    values = [float(x) for x in parts[1:]]
    rows.append([a] + values)

df = pd.DataFrame(rows, columns=["a", "E_1", "E_2", "E_3", "E_4"])

out_dir = "/mnt/data"
os.makedirs(out_dir, exist_ok=True)

file_paths = []
for i in range(1, 5):
    col = f"E_{i}"
    fig, ax = plt.subplots(figsize=(7,4.5))
    ax.plot(df["a"], df[col], marker='o')
    ax.set_yscale('log')
    ax.set_xlabel("Ruch oferowany na 1 jednostkę pojemności (a) [Erl]")
    ax.set_ylabel("Prawdopodobieństwo blokady (log scale)")
    ax.set_title(f"Prawdopodobieństwo blokady - Klasa {i} (C=250, Q=200)")
    ax.grid(True, which="both", ls="--", linewidth=0.5)
    plt.tight_layout()
    fname = f"/mnt/data/blocking_class_{i}_log.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    file_paths.append(fname)

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(df["a"], df["E_1"], marker='o', label='E_1')
ax.plot(df["a"], df["E_2"], marker='s', label='E_2')
ax.plot(df["a"], df["E_3"], marker='^', label='E_3')
ax.plot(df["a"], df["E_4"], marker='d', label='E_4')
ax.set_yscale('log')
ax.set_xlabel("Ruch oferowany na 1 jednostkę pojemności (a) [Erl]")
ax.set_ylabel("Prawdopodobieństwo blokady (log scale)")
ax.set_title("Prawdopodobieństwa blokady dla wszystkich klas (log scale)\nC=250, Q=200")
ax.grid(True, which="both", ls="--", linewidth=0.5)
ax.legend()
plt.tight_layout()
combined_fname = "/mnt/data/blocking_all_classes_log.png"
fig.savefig(combined_fname, dpi=150)
plt.close(fig)
file_paths.append(combined_fname)

table_fname = "/mnt/data/blocking_table_C250_Q200.txt"
with open(table_fname, "w") as f:
    f.write("# C = 250, Q = 200\n")
    f.write("# Columns: a E_1 E_2 E_3 E_4\n")
    for _, row in df.iterrows():
        f.write(f"{row['a']:.6f} {row['E_1']:.8e} {row['E_2']:.8e} {row['E_3']:.8e} {row['E_4']:.8e}\n")

file_paths.append(table_fname)
file_paths

