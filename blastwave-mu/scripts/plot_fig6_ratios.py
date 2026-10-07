#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE       = os.path.expanduser("~/Documents/blastwave-mu")
PHENIX_RAW = os.path.join(BASE, "input/PHENIX_spectra/AuAu/spectra.txt")
FIT        = os.path.join(BASE, "output/fits/PHENIX_AuAu/cent00/result.txt")
OUT_DIR    = os.path.join(BASE, "output/figures")
os.makedirs(OUT_DIR, exist_ok=True)

# read fit line (first non-comment line)
T_fit = bs = muB = muQ = muS = 0.0
with open(FIT) as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.split()
        if len(parts) >= 5:
            T_fit, bs, muB, muQ, muS = (float(parts[i]) for i in range(5))
            break

def parse_phenix(path, centrality=0):
    BLOCK_NAMES = ["pi+", "pi-", "K+", "K-", "p", "pbar"]
    masses = {"pi+":0.139570,"pi-":0.139570,"K+":0.493677,"K-":0.493677,
              "p":0.938272,"pbar":0.938272}
    blocks = []; cur = None
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            s = line.strip()
            if s.isdigit():
                if cur is not None: blocks.append(cur)
                cur = {"pT":[], "v":[], "e":[]}
                continue
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
    out = {}
    for name, b in zip(BLOCK_NAMES, blocks):
        pT = np.array(b["pT"]); v = np.array(b["v"]); e = np.array(b["e"])
        m = masses[name]
        out[name] = (np.sqrt(pT*pT+m*m)-m, v, e)
    return out

d = parse_phenix(PHENIX_RAW, 0)

def ratio(A,B):
    mA,vA,eA = d[A]; mB,vB,eB = d[B]
    n = min(len(mA),len(mB))
    R = vA[:n]/vB[:n]
    eR = R*np.sqrt((eA[:n]/vA[:n])**2 + (eB[:n]/vB[:n])**2)
    return mA[:n], R, eR

m1,R1,e1 = ratio("pi+","pi-")
m2,R2,e2 = ratio("K+","K-")
m3,R3,e3 = ratio("pbar","p")

mod1 = np.exp( 2*muQ/T_fit)
mod2 = np.exp( 2*(muQ+muS)/T_fit)
mod3 = np.exp(-2*(muB+muQ)/T_fit)

fig, axes = plt.subplots(1,3,figsize=(15,5))
for ax, mT, R, eR, model, ylab, col in [
    (axes[0], m1, R1, e1, mod1, r"$\pi^+/\pi^-$", "r"),
    (axes[1], m2, R2, e2, mod2, r"$K^+/K^-$",     "g"),
    (axes[2], m3, R3, e3, mod3, r"$\bar{p}/p$",   "b")]:
    ax.errorbar(mT, R, yerr=eR, fmt="o", color=col, ms=5, capsize=2, label="data")
    ax.axhline(model, color="k", ls="--", lw=2, label=f"model = {model:.3f}")
    ax.set_xlabel(r"$m_T - m_0$  [GeV]"); ax.set_ylabel(ylab)
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

fig.suptitle("PHENIX Au+Au 200 GeV, 0-5% centrality — particle/antiparticle ratios",
             fontsize=13)
plt.tight_layout(rect=[0,0,1,0.94])
plt.savefig(os.path.join(OUT_DIR, "fig6_ratios.png"), dpi=140)
print("saved fig6_ratios.png")
