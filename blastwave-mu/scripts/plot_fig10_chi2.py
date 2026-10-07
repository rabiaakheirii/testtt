#!/usr/bin/env python3
"""
Figure 10 — chi2/ndof vs N_part: the main systematic result.
Reads the existing summary.dat and Npart table.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/Documents/blastwave-mu")
FIT  = os.path.join(BASE, "output/fits/PHENIX_AuAu/summary.dat")
NP   = os.path.join(BASE, "data/Npart_PHENIX_AuAu.txt")
OUT  = os.path.join(BASE, "output/figures/fig10_chi2_vs_Npart.png")

# Load summary
d = np.loadtxt(FIT, comments="#")
cent = d[:, 0].astype(int)
chi2 = d[:, 6]

# N_part map
Nmap = {}
with open(NP) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        Nmap[int(p[0])] = float(p[1])

Npart = np.array([Nmap[c] for c in cent])
order = np.argsort(Npart)
Npart = Npart[order]
chi2  = chi2[order]

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(Npart, chi2, "o-", ms=10, lw=2.2, color="C3",
        markeredgecolor="k", markeredgewidth=0.6)

ax.axhline(1.0, color="k", ls=":", lw=1,
           label=r"Perfect fit ($\chi^2$/ndof = 1)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$N_{\rm part}$", fontsize=13)
ax.set_ylabel(r"$\chi^2 / \mathrm{ndof}$", fontsize=13)
ax.set_title("System-size dependence of the blast-wave fit quality",
             fontsize=13)
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=11)

# Annotate central & peripheral values
ax.annotate(f"{chi2[0]:.1f}",
            xy=(Npart[0], chi2[0]),
            xytext=(8, 8), textcoords="offset points",
            fontsize=11, color="darkred",
            fontweight="bold")
ax.annotate(f"{chi2[-1]:.1f}",
            xy=(Npart[-1], chi2[-1]),
            xytext=(-25, 8), textcoords="offset points",
            fontsize=11, color="darkred",
            fontweight="bold")

plt.tight_layout()
plt.savefig(OUT, dpi=140)
print(f"saved {OUT}")
