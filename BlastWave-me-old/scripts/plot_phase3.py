import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("../output/spectra/longitudinal_fit_pion.dat")
y, data, err, flow, thermal = d.T

plt.errorbar(y, data, yerr=err, fmt='o', label='data')
plt.plot(y, flow, 'r-',  label='flowing source')
plt.plot(y, thermal, 'b--', label='thermal (no flow)')
plt.xlabel('rapidity y')
plt.ylabel(r'$dn/dy$')
plt.legend()
plt.tight_layout()
plt.savefig("../output/figures/phase3_longitudinal.png", dpi=140)
print("saved")