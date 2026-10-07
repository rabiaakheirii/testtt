import numpy as np
import matplotlib.pyplot as plt
import pathlib, sys

INFILE = "../data/experimental/phenix_pion_proton_kaon.csv"
OUT    = pathlib.Path("../output/PHENIX_spectra/AuAu/figures")

species = {}
with open(INFILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        sp, x, v, e = line.split()
        species.setdefault(sp, []).append((float(x), float(v), float(e)))

fig, ax = plt.subplots(figsize=(8,6))
for sp, col in [("pi","r"), ("K","g"), ("p","b")]:
    if sp not in species:
        continue
    arr = np.array(species[sp])
    ax.errorbar(arr[:,0], arr[:,1], yerr=arr[:,2],
                fmt='o', color=col, ms=4, label=sp)

ax.set_yscale('log')
ax.set_xlabel(r'$m_T - m_0$  [GeV]')
ax.set_ylabel(r'$dN/(m_T\,dm_T)$  [GeV$^{-2}$]')
ax.set_title('PHENIX Au+Au — blast-wave input')
ax.legend()
plt.tight_layout()
plt.savefig(OUT / "phenix_spectra.png", dpi=140)
print(f"saved {OUT}/phenix_spectra.png")
