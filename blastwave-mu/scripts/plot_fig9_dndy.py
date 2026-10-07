#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.path.expanduser("~/Documents/blastwave-mu")
INFILE  = os.path.join(BASE, "output/fits/PHENIX_AuAu/cent00/dndy.dat")
OUT     = os.path.join(BASE, "output/figures/fig9_dndy.png")

species = {}
with open(INFILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        parts = line.split()
        if len(parts) < 3: continue
        sp, y, v = parts[0], float(parts[1]), float(parts[2])
        species.setdefault(sp, []).append((y, v))

fig, ax = plt.subplots(figsize=(9, 6))

style = {
    "pi+":  ("r",     "-",  r"$\pi^+$"),
    "pi-":  ("darkred","--", r"$\pi^-$"),
    "K+":   ("g",     "-",  r"$K^+$"),
    "K-":   ("darkgreen","--", r"$K^-$"),
    "p":    ("b",     "-",  r"$p$"),
    "pbar": ("navy",  "--", r"$\bar{p}$"),
}

# --- KEY CHANGE: normalize each curve to its own peak ---
for sp, (col, ls, lab) in style.items():
    if sp not in species: continue
    arr = np.array(species[sp])
    v = arr[:, 1]
    v = v / v.max()      # normalize to peak = 1
    ax.plot(arr[:, 0], v, ls, color=col, lw=2.2, label=lab)

ax.set_xlabel("rapidity  y", fontsize=13)
ax.set_ylabel("dN/dy  /  (dN/dy)$_{y=0}$", fontsize=13)
ax.set_title("Normalized rapidity distributions — PHENIX Au+Au 200 GeV, 0–5%",
             fontsize=12)
ax.grid(alpha=0.3)
ax.legend(fontsize=12, ncol=2)
ax.set_xlim(-3.5, 3.5)
plt.tight_layout()
plt.savefig(OUT, dpi=140)
print(f"saved {OUT}")
