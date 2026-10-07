import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("../output/spectra/net_proton_cumulants.dat")
sqrt_s, mu_B, cs2, C1, C2, C3, C4, Ss, ks = d.T

fig, axes = plt.subplots(2, 2, figsize=(11, 8))

ax = axes[0,0]
ax.plot(sqrt_s, cs2, 'b-o', ms=4)
ax.set_xscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ [GeV]')
ax.set_ylabel(r'$c_s^2$')
ax.set_title('Effective sound speed')
ax.grid(alpha=0.3)

ax = axes[0,1]
ax.plot(sqrt_s, mu_B, 'g-o', ms=4)
ax.set_xscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ [GeV]')
ax.set_ylabel(r'$\mu_B$ [GeV]')
ax.set_title('Baryon chemical potential')
ax.grid(alpha=0.3)

ax = axes[1,0]
ax.plot(sqrt_s, Ss, 'r-o', ms=4)
ax.axhline(0, color='k', lw=0.5)
ax.set_xscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ [GeV]')
ax.set_ylabel(r'$S\sigma = C_3/C_2$')
ax.set_title('Skewness ratio')
ax.grid(alpha=0.3)

ax = axes[1,1]
ax.plot(sqrt_s, ks, 'm-o', ms=4)
ax.axhline(0, color='k', lw=0.5)
ax.set_xscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ [GeV]')
ax.set_ylabel(r'$\kappa\sigma^2 = C_4/C_2$')
ax.set_title('Kurtosis ratio')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("../output/figures/phase8_cumulants.png", dpi=140)
print("saved phase8_cumulants.png")
