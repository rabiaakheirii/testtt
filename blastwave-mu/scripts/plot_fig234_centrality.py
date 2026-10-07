#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.path.expanduser("~/Documents/blastwave-mu")
FIT     = os.path.join(BASE, "output/fits/PHENIX_AuAu/summary.dat")
NPART   = os.path.join(BASE, "data/Npart_PHENIX_AuAu.txt")
OUT_DIR = os.path.join(BASE, "output/figures")
os.makedirs(OUT_DIR, exist_ok=True)

def load_data(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return np.array(rows)

d = load_data(FIT)
n = load_data(NPART)

print(f"[debug] fit rows = {len(d)}, npart rows = {len(n)}")

# Match by centrality index
c_fit = d[:, 0].astype(int)
c_np  = n[:, 0].astype(int)

Npart = np.zeros(len(d))
for i, ci in enumerate(c_fit):
    match = np.where(c_np == ci)[0]
    Npart[i] = n[match[0], 1] if len(match) else 0.0

T     = d[:, 1] * 1000   # MeV
bs    = d[:, 2]
muB   = d[:, 3] * 1000
muQ   = d[:, 4] * 1000
muS   = d[:, 5] * 1000
chi2  = d[:, 6]

# Figure 2
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(Npart, T, "o-", color="C0", ms=8)
ax.axhline(100, color="k", ls="--", lw=1, label="PHENIX published T ≈ 100 MeV")
ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$T_{\rm kin}$  [MeV]")
ax.set_title("Kinetic freeze-out temperature vs system size")
ax.grid(alpha=0.3); ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig2_T_vs_Npart.png"), dpi=140)
print("saved fig2_T_vs_Npart.png")

# Figure 3
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(Npart, bs, "o-", color="C3", ms=8)
ax.axhline(0.65, color="k", ls="--", lw=1, label="PHENIX published β_s ≈ 0.65")
ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\beta_s$")
ax.set_title(r"Transverse flow velocity vs system size")
ax.grid(alpha=0.3); ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig3_betas_vs_Npart.png"), dpi=140)
print("saved fig3_betas_vs_Npart.png")

# Figure 4
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(Npart, muB, "s-", color="C0", ms=7, label=r"$\mu_B$")
ax.plot(Npart, muQ, "^-", color="C2", ms=7, label=r"$\mu_Q$")
ax.plot(Npart, muS, "d-", color="C4", ms=7, label=r"$\mu_S$")
ax.axhline(0, color="k", ls=":", lw=0.8)
ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\mu$  [MeV]")
ax.set_title("Chemical potentials vs system size")
ax.grid(alpha=0.3); ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4_mu_vs_Npart.png"), dpi=140)
print("saved fig4_mu_vs_Npart.png")
