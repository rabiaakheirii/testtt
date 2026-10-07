import numpy as np
import matplotlib.pyplot as plt

for sp, m0, col in [("pi",0.1396,"r"),("K",0.4937,"g"),("p",0.9383,"b")]:
    d = np.loadtxt(f"../output/spectra/bw_{sp}.dat")
    plt.errorbar(d[:,0], d[:,1], yerr=d[:,2], fmt='o', color=col)
    plt.plot(d[:,0], d[:,3], '-', color=col, label=sp)

plt.yscale('log')
plt.xlabel(r'$m_T - m_0$ [GeV]'); plt.ylabel(r'$dn/(m_T dm_T)$')
plt.legend(); plt.tight_layout()
plt.savefig("../output/figures/phase4_blastwave.png", dpi=140)
print("saved")