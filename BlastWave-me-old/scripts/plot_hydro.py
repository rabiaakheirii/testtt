import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("../output/events/hydro_profile.dat")
tau, r, eps, T, ur = d.T

tau_vals = np.unique(tau)
cmap = plt.cm.plasma
colors = [cmap(i / (len(tau_vals)-1)) for i in range(len(tau_vals))]

# ---- energy density ----
plt.figure(figsize=(7,5))
for i, tv in enumerate(tau_vals):
    m = (tau == tv)
    plt.plot(r[m], eps[m], '-', color=colors[i], lw=1.2,
             label=r'$\tau$=%.1f' % tv if i % 5 == 0 else None)
plt.xlabel('r [fm]')
plt.ylabel(r'$\epsilon$ [GeV/fm$^3$]')
plt.title(r'Energy density vs r at different $\tau$ (Bjorken)')
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("../output/figures/phase6_energy_density.png", dpi=140)

# ---- temperature ----
plt.figure(figsize=(7,5))
for i, tv in enumerate(tau_vals):
    m = (tau == tv)
    plt.plot(r[m], T[m], '-', color=colors[i], lw=1.2,
             label=r'$\tau$=%.1f' % tv if i % 5 == 0 else None)
plt.xlabel('r [fm]')
plt.ylabel('T [GeV]')
plt.title(r'Temperature vs r at different $\tau$')
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("../output/figures/phase6_temperature.png", dpi=140)

print("saved phase6_energy_density.png and phase6_temperature.png")