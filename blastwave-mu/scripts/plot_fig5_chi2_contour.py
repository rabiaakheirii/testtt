#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.path.expanduser("~/Documents/blastwave-mu")
FIT     = os.path.join(BASE, "output/fits/PHENIX_AuAu/summary.dat")
OUT_DIR = os.path.join(BASE, "output/figures")
os.makedirs(OUT_DIR, exist_ok=True)

def load_data(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.split()
            try: rows.append([float(x) for x in parts])
            except ValueError: pass
    return np.array(rows)

d = load_data(FIT)
T  = d[:, 1] * 1000
bs = d[:, 2]
chi2 = d[:, 6]

fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(T, bs, c=chi2, cmap="viridis_r", s=150,
                edgecolor="k", linewidth=0.8,
                norm=matplotlib.colors.LogNorm())
cb = plt.colorbar(sc, label=r"$\chi^2$/ndof  (log scale)")
ax.set_xlabel("T [MeV]"); ax.set_ylabel(r"$\beta_s$")
ax.set_title(r"(T, $\beta_s$) points for 12 PHENIX centralities, colored by $\chi^2$/ndof")
ax.grid(alpha=0.3)
# annotate centrality indices
for i, ci in enumerate(d[:, 0].astype(int)):
    ax.annotate(f"{ci}", (T[i], bs[i]),
                xytext=(4, 4), textcoords="offset points", fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig5_chi2_contour.png"), dpi=140)
print("saved fig5_chi2_contour.png")
