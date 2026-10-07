#!/usr/bin/env python3
"""Compare standard blast-wave results with published PHENIX values."""
import os
import numpy as np

BASE    = os.path.expanduser("~/Documents/blastwave-mu")
SUMMARY = os.path.join(BASE, "output/fits/PHENIX_AuAu_standard/summary.dat")
NPART   = os.path.join(BASE, "data/Npart_PHENIX_AuAu.txt")

# Published PHENIX values (PRC 69, 034909)
# for different centralities, roughly:
published = {
    # N_part: (T_kin_MeV, beta_s)
    351: (100, 0.62),   # 0-5%
    114: (105, 0.62),   # 30-40%
    22:  (115, 0.65),   # 70-80%
    7:   (120, 0.70),   # 90-100%
}

d = np.loadtxt(SUMMARY, comments="#")
cent  = d[:, 0].astype(int)
T     = d[:, 1] * 1000
bs    = d[:, 2]
muB   = d[:, 3] * 1000
chi2  = d[:, 6]

# N_part map
Nmap = {}
with open(NPART) as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        p = line.split()
        Nmap[int(p[0])] = float(p[1])

print()
print("="*80)
print(" STANDARD BLAST-WAVE RESULTS (μ_Q = μ_S = 0)")
print("="*80)
print()
print(f"{'cent':>4}  {'N_part':>7}  {'T [MeV]':>8}  {'β_s':>6}  {'μ_B [MeV]':>10}  {'χ²/ndof':>8}  {'vs PHENIX':>10}")
print("-" * 80)

for i in range(len(cent)):
    c = cent[i]
    N = Nmap.get(c, 0)
    # Find nearest published value
    closest = min(published.keys(), key=lambda x: abs(x - N))
    T_pub, bs_pub = published[closest]
    dT = T[i] - T_pub
    dbs = bs[i] - bs_pub
    match = "✅" if abs(dT) < 20 and abs(dbs) < 0.10 else \
            "⚠️" if abs(dT) < 40 and abs(dbs) < 0.20 else "❌"
    print(f"{c:>4}  {N:>7.0f}  {T[i]:>8.1f}  {bs[i]:>6.3f}  {muB[i]:>10.1f}  {chi2[i]:>8.1f}  {match:>10}")

print()
print("="*80)
print(" PASS/FAIL CRITERIA")
print("="*80)
print("  Peripheral (N_part < 30): β_s should be 0.65-0.75, T 100-140 MeV")
print("  Central (N_part > 200):   β_s should be 0.55-0.65, T 90-120 MeV")
print()
