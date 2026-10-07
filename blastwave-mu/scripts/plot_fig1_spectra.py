#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.path.expanduser("~/Documents/blastwave-mu")
FIT_DIR = os.path.join(BASE, "output/fits/PHENIX_AuAu/cent00")
OUT_DIR = os.path.join(BASE, "output/figures")
os.makedirs(OUT_DIR, exist_ok=True)

SPECIES = [
    ("pi+",  r"$\pi^+$", "r",  "o"),
    ("pi-",  r"$\pi^-$", "r",  "s"),
    ("K+",   r"$K^+$",   "g",  "o"),
    ("K-",   r"$K^-$",   "g",  "s"),
    ("p",    r"$p$",     "b",  "o"),
    ("pbar", r"$\bar{p}$","b", "s"),
]

fig, axes = plt.subplots(3, 2, figsize=(12, 12))
fig.suptitle("PHENIX Au+Au 200 GeV, 0-5% centrality\n"
             r"Blast-wave fit with $\mu_B$, $\mu_Q$, $\mu_S$ and 34-resonance feed-down",
             fontsize=13)

for idx, (name, label, color, marker) in enumerate(SPECIES):
    ax = axes[idx // 2, idx % 2]
    fname = os.path.join(FIT_DIR, f"fit_{name}.dat")
    if not os.path.exists(fname):
        ax.text(0.5, 0.5, "missing", ha="center", va="center",
                transform=ax.transAxes); continue
    d = np.loadtxt(fname, comments="#")
    ax.errorbar(d[:, 0], d[:, 1], yerr=d[:, 2], fmt=marker, color=color,
                ms=5, capsize=2, label="data")
    ax.plot(d[:, 0], d[:, 3], "k-", lw=2, label="model")
    ax.set_yscale("log")
    ax.set_xlabel(r"$m_T - m_0$  [GeV]")
    ax.set_ylabel(r"$dN/(m_T\,dm_T)$  [GeV$^{-2}$]")
    ax.set_title(label)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(OUT_DIR, "fig1_spectra.png"), dpi=140)
print("saved fig1_spectra.png")
