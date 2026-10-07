#!/usr/bin/env python3
"""
Generate all 7 figures from the blastwave-mu pipeline.

Data sources:
  output/fits/PHENIX_AuAu/cent00/result.txt         (Phase 7 best fit)
  output/fits/PHENIX_AuAu/cent00/fit_*.dat          (Phase 7 best-fit curves)
  output/fits/PHENIX_AuAu/cent00/chi2_surface.dat   (Phase 10)
  output/fits/PHENIX_AuAu/cent00/posterior_stats.dat(Phase 12)
  output/fits/PHENIX_AuAu/params_vs_ptbin.dat       (Phase 11)
  output/fits/PHENIX_AuAu/summary.dat               (Phase 9)
  data/Npart_PHENIX_AuAu.txt                        (N_part)
  input/PHENIX_spectra/AuAu/spectra.txt             (raw data, Fig 8)

Outputs:
  output/figures/fig1_globalfit.png
  output/figures/fig2_contour.png
  output/figures/fig3_T_vs_ptbin.png
  output/figures/fig4_ut_vs_ptbin.png
  output/figures/fig5_T_vs_Npart.png
  output/figures/fig6_betas_vs_Npart.png
  output/figures/fig7_mu_vs_Npart.png
  output/figures/fig8_ratios.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE    = os.path.expanduser("~/Documents/blastwave-mu")
FITDIR  = os.path.join(BASE, "output/fits/PHENIX_AuAu")
CENT00  = os.path.join(FITDIR, "cent00")
OUT     = os.path.join(BASE, "output/figures")
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------
# Utility
# ------------------------------------------------------------------
def load_ascii(path, **kw):
    if not os.path.exists(path):
        print(f"  [skip] {path} not found")
        return None
    return np.loadtxt(path, comments="#", **kw)

def load_2col(path):
    d = load_ascii(path)
    if d is None: return None
    return d[:, 0], d[:, 1]

# ==================================================================
# Figure 1 — Global fit: 6 panels of π, K, p spectra
# ==================================================================
print("[1/8] Figure 1: global fit")

SPECIES = [
    ("pi+",  r"$\pi^+$",  "red",    "o"),
    ("pi-",  r"$\pi^-$",  "red",    "s"),
    ("K+",   r"$K^+$",    "green",  "o"),
    ("K-",   r"$K^-$",    "green",  "s"),
    ("p",    r"$p$",      "blue",   "o"),
    ("pbar", r"$\bar{p}$","blue",   "s"),
]

fig, axes = plt.subplots(3, 2, figsize=(12, 12))
fig.suptitle("PHENIX Au+Au 200 GeV, 0–5% centrality\n"
             r"Blast-wave with $\mu_B,\mu_Q,\mu_S$ + 34-resonance feed-down",
             fontsize=13, y=0.995)

for idx, (name, label, color, marker) in enumerate(SPECIES):
    ax = axes[idx // 2, idx % 2]
    fname = os.path.join(CENT00, f"fit_{name}.dat")
    if not os.path.exists(fname):
        ax.text(0.5, 0.5, "missing", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(label); continue

    d = np.loadtxt(fname, comments="#")
    mT, val, err, model = d[:, 0], d[:, 1], d[:, 2], d[:, 3]

    ax.errorbar(mT, val, yerr=err, fmt=marker, color=color,
                ms=5, capsize=2, label="data", zorder=2)
    ax.plot(mT, model, "k-", lw=2, label="model", zorder=3)

    ax.set_yscale("log")
    ax.set_xlabel(r"$m_T - m_0$  [GeV]")
    ax.set_ylabel(r"$dN/(m_T\,dm_T)$  [GeV$^{-2}$]")
    ax.set_title(label, fontsize=14)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper right")

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(os.path.join(OUT, "fig1_globalfit.png"), dpi=140)
plt.close()
print("    saved fig1_globalfit.png")

# ==================================================================
# Figure 2 — (T, β_s) χ² contour
# ==================================================================
print("[2/8] Figure 2: χ² contour")

d = load_ascii(os.path.join(CENT00, "chi2_surface.dat"))
if d is not None:
    T   = d[:, 0] * 1000    # MeV
    bs  = d[:, 1]
    chi2 = d[:, 2]

    Tn = np.unique(T)
    bn = np.unique(bs)
    Z  = chi2.reshape(len(bn), len(Tn))
    dchi2 = Z - Z.min()

    fig, ax = plt.subplots(figsize=(8, 6))
    cf = ax.contourf(Tn, bn, np.log10(dchi2 + 1e-3),
                     levels=25, cmap="viridis")
    plt.colorbar(cf, label=r"$\log_{10}(\Delta\chi^2)$")

    # 1σ, 2σ, 3σ for 2 parameters
    cs = ax.contour(Tn, bn, dchi2,
                    levels=[2.30, 6.18, 11.83],
                    colors=["white", "yellow", "red"],
                    linewidths=1.8)
    fmt = {2.30: r"$1\sigma$", 6.18: r"$2\sigma$", 11.83: r"$3\sigma$"}
    ax.clabel(cs, fmt=fmt, fontsize=11, inline=True)

    # Best-fit marker
    imin = np.unravel_index(np.argmin(Z), Z.shape)
    ax.plot(Tn[imin[1]], bn[imin[0]], "w*", ms=20,
            mec="k", mew=1, label="best fit")

    ax.set_xlabel("T [MeV]", fontsize=12)
    ax.set_ylabel(r"$\beta_s$", fontsize=12)
    ax.set_title(r"$\chi^2$ contour in ($T$, $\beta_s$) — centrality 0–5%",
                 fontsize=13)
    ax.legend(fontsize=11, loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig2_contour.png"), dpi=140)
    plt.close()
    print("    saved fig2_contour.png")

# ==================================================================
# Figure 3 — T vs pT window
# ==================================================================
print("[3/8] Figure 3: T vs pT window")

d = load_ascii(os.path.join(FITDIR, "params_vs_ptbin.dat"))
if d is not None:
    pT_min = d[:, 1]
    pT_max = d[:, 2]
    T_val  = d[:, 3] * 1000
    bs_val = d[:, 4]
    T_avg  = 0.5 * (pT_min + pT_max)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(T_avg, T_val, "o-", ms=10, lw=2, color="C0")
    ax.axhline(100, color="k", ls="--", lw=1,
               label="Published PHENIX: T = 100 MeV")
    ax.set_xlabel(r"$\langle p_T \rangle$ window center  [GeV]", fontsize=12)
    ax.set_ylabel(r"$T_{\rm kin}$  [MeV]", fontsize=12)
    ax.set_title(r"Extracted $T_{\rm kin}$ vs $p_T$ window", fontsize=13)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig3_T_vs_ptbin.png"), dpi=140)
    plt.close()
    print("    saved fig3_T_vs_ptbin.png")

# ==================================================================
# Figure 4 — β_s vs pT window
# ==================================================================
print("[4/8] Figure 4: β_s vs pT window")

if d is not None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(T_avg, bs_val, "s-", ms=10, lw=2, color="C3")
    ax.axhline(0.62, color="k", ls="--", lw=1,
               label=r"Published PHENIX: $\beta_s$ = 0.62")
    ax.set_xlabel(r"$\langle p_T \rangle$ window center  [GeV]", fontsize=12)
    ax.set_ylabel(r"$\beta_s$", fontsize=12)
    ax.set_title(r"Extracted $\beta_s$ vs $p_T$ window", fontsize=13)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig4_ut_vs_ptbin.png"), dpi=140)
    plt.close()
    print("    saved fig4_ut_vs_ptbin.png")

# ==================================================================
# Figure 5 — T vs N_part
# ==================================================================
print("[5/8] Figure 5: T vs N_part")

d_sum = load_ascii(os.path.join(FITDIR, "summary.dat"))
d_np  = load_ascii(os.path.join(BASE, "data/Npart_PHENIX_AuAu.txt"))

if d_sum is not None and d_np is not None:
    cent  = d_sum[:, 0].astype(int)
    T_c   = d_sum[:, 1] * 1000
    bs_c  = d_sum[:, 2]
    muB_c = d_sum[:, 3] * 1000
    muQ_c = d_sum[:, 4] * 1000
    muS_c = d_sum[:, 5] * 1000
    chi2  = d_sum[:, 6]

    Nmap = dict(zip(d_np[:, 0].astype(int), d_np[:, 1]))
    Npart = np.array([Nmap[c] for c in cent])

    # sort by Npart
    order = np.argsort(Npart)
    Npart = Npart[order]
    T_c   = T_c[order]
    bs_c  = bs_c[order]
    muB_c = muB_c[order]
    muQ_c = muQ_c[order]
    muS_c = muS_c[order]
    chi2  = chi2[order]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(Npart, T_c, "o-", ms=10, lw=2, color="C0")
    ax.axhline(100, color="k", ls="--", lw=1,
               label="Published PHENIX: T = 100 MeV")
    ax.set_xscale("log")
    ax.set_xlabel(r"$N_{\rm part}$", fontsize=12)
    ax.set_ylabel(r"$T_{\rm kin}$  [MeV]", fontsize=12)
    ax.set_title(r"$T_{\rm kin}$ vs system size (12 centralities)", fontsize=13)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig5_T_vs_Npart.png"), dpi=140)
    plt.close()
    print("    saved fig5_T_vs_Npart.png")

# ==================================================================
# Figure 6 — β_s vs N_part
# ==================================================================
print("[6/8] Figure 6: β_s vs N_part")

if d_sum is not None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(Npart, bs_c, "s-", ms=10, lw=2, color="C3")
    ax.axhline(0.62, color="k", ls="--", lw=1,
               label=r"Published PHENIX: $\beta_s$ = 0.62")
    ax.set_xscale("log")
    ax.set_xlabel(r"$N_{\rm part}$", fontsize=12)
    ax.set_ylabel(r"$\beta_s$", fontsize=12)
    ax.set_title(r"$\beta_s$ vs system size", fontsize=13)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig6_betas_vs_Npart.png"), dpi=140)
    plt.close()
    print("    saved fig6_betas_vs_Npart.png")

# ==================================================================
# Figure 7 — μ_B, μ_Q, μ_S vs N_part
# ==================================================================
print("[7/8] Figure 7: chemical potentials vs N_part")

if d_sum is not None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(Npart, muB_c, "s-", ms=9, lw=2, color="C0",
            label=r"$\mu_B$")
    ax.plot(Npart, muQ_c, "^-", ms=9, lw=2, color="C2",
            label=r"$\mu_Q$")
    ax.plot(Npart, muS_c, "d-", ms=9, lw=2, color="C4",
            label=r"$\mu_S$")
    ax.axhline(0, color="k", ls=":", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"$N_{\rm part}$", fontsize=12)
    ax.set_ylabel(r"$\mu$  [MeV]", fontsize=12)
    ax.set_title(r"Chemical potentials vs system size", fontsize=13)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig7_mu_vs_Npart.png"), dpi=140)
    plt.close()
    print("    saved fig7_mu_vs_Npart.png")

# ==================================================================
# Figure 8 — particle/antiparticle ratios
# ==================================================================
print("[8/8] Figure 8: particle/antiparticle ratios")

RAW = os.path.join(BASE, "input/PHENIX_spectra/AuAu/spectra.txt")
RES = os.path.join(CENT00, "result.txt")

def parse_phenix_raw(path, centrality=0):
    BLOCKS = ["pi+", "pi-", "K+", "K-", "p", "pbar"]
    MASSES = {"pi+":0.139570,"pi-":0.139570,"K+":0.493677,
              "K-":0.493677,"p":0.938272,"pbar":0.938272}
    blocks = []; cur = None
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            s = line.strip()
            if s.isdigit():
                if cur is not None: blocks.append(cur)
                cur = {"pT":[], "v":[], "e":[]}; continue
            if cur is None: continue
            parts = line.split()
            try: pT = float(parts[0])
            except ValueError: continue
            pairs = []
            for i in range(1, len(parts)-1, 2):
                try: pairs.append((float(parts[i]), abs(float(parts[i+1]))))
                except ValueError: break
            if centrality < len(pairs):
                v,e = pairs[centrality]
                cur["pT"].append(pT); cur["v"].append(v); cur["e"].append(e)
    if cur is not None: blocks.append(cur)
    result = {}
    for name, b in zip(BLOCKS, blocks):
        pT = np.array(b["pT"]); v = np.array(b["v"]); e = np.array(b["e"])
        m = MASSES[name]
        mT = np.sqrt(pT*pT + m*m) - m
        result[name] = (mT, v, e)
    return result

if os.path.exists(RAW):
    data = parse_phenix_raw(RAW, 0)

    # read fit params
    T_fit = bs_fit = muB = muQ = muS = 0.0
    with open(RES) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) >= 5:
                T_fit, bs_fit, muB, muQ, muS = (float(parts[i]) for i in range(5))
                break

    def ratio(A, B):
        mA, vA, eA = data[A]; mB, vB, eB = data[B]
        n = min(len(mA), len(mB))
        R = vA[:n] / vB[:n]
        eR = R * np.sqrt((eA[:n]/vA[:n])**2 + (eB[:n]/vB[:n])**2)
        return mA[:n], R, eR

    mT1, R1, e1 = ratio("pi+", "pi-")
    mT2, R2, e2 = ratio("K+",  "K-")
    mT3, R3, e3 = ratio("pbar","p")

    mod1 = np.exp( 2.0 * muQ / T_fit)
    mod2 = np.exp( 2.0 * (muQ + muS) / T_fit)
    mod3 = np.exp(-2.0 * (muB + muQ) / T_fit)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, mT, R, eR, model, ylab, col in [
        (axes[0], mT1, R1, e1, mod1, r"$\pi^+/\pi^-$", "r"),
        (axes[1], mT2, R2, e2, mod2, r"$K^+/K^-$",     "g"),
        (axes[2], mT3, R3, e3, mod3, r"$\bar{p}/p$",   "b")]:
        ax.errorbar(mT, R, yerr=eR, fmt="o", color=col, ms=6,
                    capsize=2, label="data")
        ax.axhline(model, color="k", ls="--", lw=2,
                   label=f"model = {model:.3f}")
        ax.set_xlabel(r"$m_T - m_0$  [GeV]", fontsize=12)
        ax.set_ylabel(ylab, fontsize=13)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=11)

    fig.suptitle("PHENIX Au+Au 200 GeV, 0–5% centrality — "
                 "particle/antiparticle ratios", fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(OUT, "fig8_ratios.png"), dpi=140)
    plt.close()
    print("    saved fig8_ratios.png")

print("\nAll figures written to:", OUT)
