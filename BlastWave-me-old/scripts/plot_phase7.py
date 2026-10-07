import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("../output/spectra/cooper_frye_pion.dat")
mT, dN = d.T

plt.figure(figsize=(7,5))
plt.plot(mT, dN, 'b-', lw=1.5, label='Cooper-Frye (Bjorken)')
plt.yscale('log')
plt.xlabel(r'$m_T$ [GeV]')
plt.ylabel(r'$dN/(m_T\,dm_T)$')
plt.title(r'Pion $m_T$ spectrum from Cooper-Frye')
plt.legend()
plt.tight_layout()
plt.savefig("../output/figures/phase7_cooper_frye.png", dpi=140)
print("saved phase7_cooper_frye.png")
