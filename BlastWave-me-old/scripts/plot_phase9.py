import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# --- (0,0) pi, K, p blast-wave fits ---
ax = axes[0, 0]
for sp, col in [("pi", "r"), ("K", "g"), ("p", "b")]:
    try:
        d = np.loadtxt(f"../output/spectra/bw_{sp}.dat")
        ax.errorbar(d[:, 0], d[:, 1], yerr=d[:, 2],
                    fmt='o', color=col, ms=3)
        ax.plot(d[:, 0], d[:, 3], '-', color=col, label=sp)
    except OSError:
        pass
ax.set_yscale('log')
ax.set_xlabel('mT - m0  [GeV]')
ax.set_ylabel('dN / (mT dmT)')
ax.set_title('Blast-wave fits (Phase 4)')
ax.legend()

# --- (0,1) thermal data vs Cooper-Frye ---
ax = axes[0, 1]
try:
    d1 = np.loadtxt("../output/spectra/thermal_fit_pion.dat")
    ax.errorbar(d1[:, 0], d1[:, 1], yerr=d1[:, 2],
                fmt='o', ms=3, label='data')
    ax.plot(d1[:, 0], d1[:, 3], 'r-', label='thermal fit')
except OSError:
    pass
try:
    d2 = np.loadtxt("../output/spectra/cooper_frye_pion.dat")
    scale = d1[:, 1].max() / d2[:, 1].max()
    ax.plot(d2[:, 0] - 0.13957, d2[:, 1] * scale,
            'b--', label='Cooper-Frye')
except OSError:
    pass
ax.set_yscale('log')
ax.set_xlabel('mT - m0  [GeV]')
ax.set_ylabel('dN / (mT dmT)')
ax.set_title('Thermal vs Cooper-Frye')
ax.legend()

# --- (1,0) cs2 vs sqrt(s) ---
ax = axes[1, 0]
try:
    d = np.loadtxt("../output/spectra/net_proton_cumulants.dat")
    sqrt_s, cs2 = d[:, 0], d[:, 2]
    ax.plot(sqrt_s, cs2, 'b-o', ms=4)
    ax.set_xscale('log')
    ax.set_xlabel('sqrt(s_NN)  [GeV]')
    ax.set_ylabel('cs2')
    ax.set_title('Sound speed vs collision energy (Phase 8)')
    ax.grid(alpha=0.3)
except OSError:
    pass

# --- (1,1) energy density vs r ---
ax = axes[1, 1]
try:
    h = np.loadtxt("../output/events/hydro_profile.dat")
    tau, r, eps = h[:, 0], h[:, 1], h[:, 2]
    tau_vals = np.unique(tau)
    step = max(1, len(tau_vals) // 8)
    for tv in tau_vals[::step]:
        m = (tau == tv)
        ax.plot(r[m], eps[m], '-', lw=1.2, label=f"tau={tv:.1f}")
    ax.set_xlabel('r  [fm]')
    ax.set_ylabel('epsilon  [GeV/fm^3]')
    ax.set_title('Hydro: energy density (Phase 6)')
    ax.legend(fontsize=7)
except OSError:
    pass

plt.tight_layout()
plt.savefig("../output/figures/phase9_summary.png", dpi=140)
print("saved ../output/figures/phase9_summary.png")
