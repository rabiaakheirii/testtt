#!/usr/bin/env python3
"""
NICA/MPD 9.2 GeV analysis plots.

Produces:
  nica_fig1_spectra.png    -- spectra + fits for 3 centralities
  nica_fig2_params_vs_Npart.png -- T, beta_s, mu_B, mu_Q, mu_S vs N_part
  nica_fig3_ratios.png     -- pi+/pi-, K+/K-, pbar/p vs mT
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/Documents/blastwave-mu")
FITD = os.path.join(BASE, "output/fits/NICA_MPD")
OUT  = os.path.join(BASE, "output/figures")
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------
# Figure 1: spectra + fits for 3 centralities
# ------------------------------------------------------------------
SPECIES = [
    ("pi+",  r"$\pi^+$",  "red",    "o"),
    ("K+",   r"$K^+$",    "green",  "s"),
    ("p",    r"$p$",      "blue",   "D"),
    ("d",    r"$d$",      "orange", "^"),
]

cent_labels = ["0-20%", "20-40%", "40-80%"]

fig, axes = plt.subplots(len(SPECIES), 3, figsize=(13, 12))
fig.suptitle("NICA/MPD Au+Au 9.2 GeV\nBlast-wave fit with $\\mu_B, \\mu_Q, \\mu_S$",
             fontsize=13, y=0.995)

for row, (name, label, color, marker) in enumerate(SPECIES):
    for col in range(3):
        ax = axes[row, col]
        fname = os.path.join(FITD, f"cent{col}", f"fit_{name}.dat")
        if not os.path.exists(fname):
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes)
            if row == 0: ax.set_title(cent_labels[col])
            continue
        d = np.loadtxt(fname, comments="#")
        if d.ndim < 2 or len(d) == 0:
            continue
        mT, val, err, model = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
        ax.errorbar(mT, val, yerr=err, fmt=marker, color=color,
                    ms=5, capsize=2, label="data")
        ax.plot(mT, model, "k-", lw=2, label="model")
        ax.set_yscale("log")
        ax.grid(alpha=0.3)
        if row == 0: ax.set_title(cent_labels[col], fontsize=12)
        if col == 0: ax.set_ylabel(f"{label}\n$dN/(m_T\\,dm_T)$", fontsize=11)
        if row == len(SPECIES) - 1: ax.set_xlabel(r"$m_T - m_0$ [GeV]")
        if row == 0 and col == 0:
            ax.legend(fontsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(OUT, "nica_fig1_spectra.png"), dpi=140)
plt.close()
print("saved nica_fig1_spectra.png")

# ------------------------------------------------------------------
# Figure 2: parameters vs N_part
# ------------------------------------------------------------------
summary = os.path.join(FITD, "summary.dat")
if os.path.exists(summary):
    d = np.loadtxt(summary, comments="#")
    if d.ndim == 1:
        d = d.reshape(1, -1)
    bin_idx, Npart = d[:, 0].astype(int), d[:, 1]
    T      = d[:, 2] * 1000
    bs     = d[:, 3]
    muB    = d[:, 4] * 1000
    muQ    = d[:, 5] * 1000
    muS    = d[:, 6] * 1000
    chi2   = d[:, 7]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("NICA/MPD 9.2 GeV — extracted parameters", fontsize=13)

    ax = axes[0, 0]
    ax.plot(Npart, T, "o-", ms=10, lw=2, color="C0")
    ax.axhline(130, color="k", ls="--", lw=1, label="STAR BES 9.2 GeV: ~130 MeV")
    ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel("T [MeV]")
    ax.set_title("Kinetic freeze-out T"); ax.grid(alpha=0.3); ax.legend()

    ax = axes[0, 1]
    ax.plot(Npart, bs, "s-", ms=10, lw=2, color="C3")
    ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\beta_s$")
    ax.set_title("Surface flow velocity"); ax.grid(alpha=0.3)

    ax = axes[0, 2]
    ax.plot(Npart, muB, "D-", ms=10, lw=2, color="C2")
    ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\mu_B$ [MeV]")
    ax.set_title(r"Baryon chemical potential")
    ax.grid(alpha=0.3)
    ax.axhline(420, color="k", ls="--", lw=1, label="Expected ~420 MeV")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(Npart, muQ, "^-", ms=10, lw=2, color="C4")
    ax.axhline(0, color="k", ls=":", lw=0.8)
    ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\mu_Q$ [MeV]")
    ax.set_title(r"Electric charge potential"); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(Npart, muS, "d-", ms=10, lw=2, color="C5")
    ax.axhline(0, color="k", ls=":", lw=0.8)
    ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\mu_S$ [MeV]")
    ax.set_title(r"Strangeness potential"); ax.grid(alpha=0.3)

    ax = axes[1, 2]
    ax.plot(Npart, chi2, "o-", ms=10, lw=2, color="C1")
    ax.set_yscale("log")
    ax.set_xlabel(r"$N_{\rm part}$"); ax.set_ylabel(r"$\chi^2/{\rm ndof}$")
    ax.set_title("Fit quality"); ax.grid(alpha=0.3, which="both")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(OUT, "nica_fig2_params_vs_Npart.png"), dpi=140)
    plt.close()
    print("saved nica_fig2_params_vs_Npart.png")

    # Print summary
    print("\n=== NICA Summary ===")
    print(f"{'bin':>4} {'N_part':>7} {'T[MeV]':>8} {'β_s':>6} "
          f"{'μ_B[MeV]':>10} {'μ_Q[MeV]':>10} {'μ_S[MeV]':>10} {'χ²/ndof':>10}")
    for i in range(len(bin_idx)):
        print(f"{bin_idx[i]:>4} {Npart[i]:>7.0f} {T[i]:>8.1f} {bs[i]:>6.3f} "
              f"{muB[i]:>10.1f} {muQ[i]:>10.1f} {muS[i]:>10.1f} {chi2[i]:>10.2f}")

# ------------------------------------------------------------------
# Figure 3: particle ratios
# ------------------------------------------------------------------
def parse_mpd(path, centrality=0):
    """Parse MPD file returning {species: (mT, value, err)}."""
    data = {}
    current_species = None
    current_cent = -1
    m = {}
    with open(path) as f:
        for line in f:
            if line.startswith("# centrality"):
                current_cent += 1
                continue
            if line.startswith("# particle:"):
                current_species = line.split(":")[1].strip()
                continue
            if line.startswith("#") or not line.strip():
                continue
            if current_cent != centrality: continue
            parts = line.split()
            if len(parts) < 3: continue
            try:
                pT, val, err = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            if current_species not in data:
                data[current_species] = {"pT": [], "val": [], "err": []}
            data[current_species]["pT"].append(pT)
            data[current_species]["val"].append(val)
            data[current_species]["err"].append(err)
    # Convert to numpy
    out = {}
    for sp, d in data.items():
        # mass map
        mass = {"pi+":0.1396,"pi-":0.1396,"K+":0.4937,"K-":0.4937,
                "p":0.9383,"pbar":0.9383,"d":1.8756,"t":2.8089,
                "He3":2.8089,"He4":3.7274}.get(sp, 0.1396)
        pT = np.array(d["pT"])
        mT = np.sqrt(pT*pT + mass*mass) - mass
        out[sp] = (mT, np.array(d["val"]), np.array(d["err"]))
    return out

MPD_RAW = os.path.join(BASE, "input/MPD_sim/AuAu/mpd_spectra.txt")
if os.path.exists(MPD_RAW):
    data = parse_mpd(MPD_RAW, 0)

    def ratio(A, B):
        mA, vA, eA = data[A]
        mB, vB, eB = data[B]
        n = min(len(mA), len(mB))
        R = vA[:n] / vB[:n]
        eR = R * np.sqrt((eA[:n]/vA[:n])**2 + (eB[:n]/vB[:n])**2)
        return mA[:n], R, eR

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (A, B), color, ylabel in [
        (axes[0], ("pi+","pi-"), "r", r"$\pi^+/\pi^-$"),
        (axes[1], ("K+","K-"),   "g", r"$K^+/K^-$"),
        (axes[2], ("pbar","p"),  "b", r"$\bar{p}/p$")]:
        if A in data and B in data:
            mT, R, eR = ratio(A, B)
            ax.errorbar(mT, R, yerr=eR, fmt="o", color=color,
                        ms=6, capsize=2)
            ax.set_xlabel(r"$m_T - m_0$ [GeV]"); ax.set_ylabel(ylabel)
            ax.set_title(ylabel); ax.grid(alpha=0.3)

    fig.suptitle("NICA/MPD 9.2 GeV, centrality 0-20% — particle/antiparticle ratios",
                 fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(OUT, "nica_fig3_ratios.png"), dpi=140)
    plt.close()
    print("saved nica_fig3_ratios.png")

print("\nAll NICA figures written to:", OUT)
