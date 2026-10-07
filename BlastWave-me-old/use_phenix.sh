#!/usr/bin/env bash
set -e
cd ~/Documents/qcd-collision-sim

# -------------------- paths --------------------
INPUT_DIR="$HOME/Documents/qcd-collision-sim/input/PHENIX_spectra/AuAu"
OUTPUT_DIR="$HOME/Documents/qcd-collision-sim/output/PHENIX_spectra/AuAu"

CENTRALITY_COL="${1:-0}"        # 0 = most central, 11 = most peripheral

echo "=================================================="
echo " PHENIX Au+Au -> blast-wave pipeline"
echo " input  : $INPUT_DIR"
echo " output : $OUTPUT_DIR"
echo " centrality column: $CENTRALITY_COL"
echo "=================================================="

# -------------------- sanity checks --------------------
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: input directory does not exist: $INPUT_DIR"
    exit 1
fi

# find the data file (any regular file inside INPUT_DIR)
PHENIX_SRC=""
for f in "$INPUT_DIR"/*; do
    if [ -f "$f" ]; then
        PHENIX_SRC="$f"
        break
    fi
done

if [ -z "$PHENIX_SRC" ]; then
    echo "ERROR: no data file found in $INPUT_DIR"
    exit 1
fi

echo "data file: $PHENIX_SRC"

# -------------------- create output tree --------------------
mkdir -p "$OUTPUT_DIR"/{spectra,figures,logs}
mkdir -p data/experimental

# -------------------- 1. parse PHENIX file --------------------
PHENIX_SRC="$PHENIX_SRC" CENTRALITY_COL="$CENTRALITY_COL" OUTPUT_DIR="$OUTPUT_DIR" python3 - <<'PYEOF'
import os, math, pathlib

src       = os.environ["PHENIX_SRC"]
col       = int(os.environ["CENTRALITY_COL"])
outdir    = pathlib.Path(os.environ["OUTPUT_DIR"])

MASSES = {"pi": 0.13957, "K": 0.49367, "p": 0.93827}
BLOCK_SPECIES = {0: "pi", 1: "pi", 2: "K", 3: "K", 4: "p", 5: "p"}

def parse(path):
    blocks, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line.strip():
                continue
            parts = line.split()
            # a line with a single integer => new block header
            if len(parts) == 1 and parts[0].isdigit():
                if cur is not None:
                    blocks.append(cur)
                cur = {"n": int(parts[0]), "rows": []}
                continue
            if cur is None:
                continue
            pT = float(parts[0])
            pairs = []
            for i in range(1, len(parts) - 1, 2):
                try:
                    pairs.append((float(parts[i]), float(parts[i+1])))
                except ValueError:
                    break
            cur["rows"].append((pT, pairs))
    if cur is not None:
        blocks.append(cur)
    return blocks

blocks = parse(src)
print(f"[parse] {len(blocks)} blocks found in {src}")
for i, b in enumerate(blocks):
    nc = len(b["rows"][0][1]) if b["rows"] else 0
    print(f"   block {i+1}: header={b['n']:3d}  rows={len(b['rows']):3d}  centrality bins={nc}")

# write output file
outfile = pathlib.Path("data/experimental/phenix_pion_proton_kaon.csv")
lines = [f"# PHENIX Au+Au — centrality column {col}",
         "# species  mT-m0[GeV]  dN/(mT dmT)  err"]

used_species = set()
n_pts = 0
for bi, block in enumerate(blocks):
    sp = BLOCK_SPECIES.get(bi)
    if sp is None or sp in used_species:
        continue
    used_species.add(sp)
    m = MASSES[sp]
    for pT, pairs in block["rows"]:
        if col >= len(pairs):
            continue
        val, err = pairs[col]
        if val <= 0:
            continue
        mT = math.sqrt(pT*pT + m*m)
        lines.append(f"{sp}  {mT-m:.6e}  {val:.6e}  {err:.6e}")
        n_pts += 1

outfile.write_text("\n".join(lines) + "\n")
print(f"[convert] wrote {outfile}  ({n_pts} points, species: {sorted(used_species)})")

# also copy raw + parsed to output dir
import shutil
shutil.copy(src, outdir / "spectra" / "phenix_raw.txt")
shutil.copy(outfile, outdir / "spectra" / "phenix_parsed.csv")
print(f"[io] copied raw + parsed files to {outdir}/spectra/")
PYEOF

# -------------------- 2. patch main.cpp --------------------
if [ ! -f data/experimental/na35_pion_proton_kaon.csv.bak ]; then
    cp data/experimental/na35_pion_proton_kaon.csv \
       data/experimental/na35_pion_proton_kaon.csv.bak 2>/dev/null || true
fi

python3 - <<'PYEOF'
import pathlib
p = pathlib.Path("src/main.cpp")
t = p.read_text()
t = t.replace("na35_pion_proton_kaon.csv", "phenix_pion_proton_kaon.csv")
p.write_text(t)
print("[patch] main.cpp now uses phenix_pion_proton_kaon.csv")
PYEOF

# -------------------- 3. preview --------------------
echo ""
echo "---- preview ----"
head -6 data/experimental/phenix_pion_proton_kaon.csv
echo "..."
grep -c '^[piKp]' data/experimental/phenix_pion_proton_kaon.csv

# -------------------- 4. build & run --------------------
echo ""
echo "---- building ----"
cd build && cmake .. >/dev/null && make && cd ..

echo ""
echo "---- running ----"
./build/sim 2>&1 | tee "$OUTPUT_DIR/logs/run.log" | grep -A 4 "Blast-wave fit"

# -------------------- 5. plot --------------------
echo ""
echo "---- plotting ----"
cat > scripts/plot_phenix_bw.py <<'EOF'
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
EOF
cd scripts && python3 plot_phenix_bw.py && cd ..

# -------------------- 6. copy blast-wave fits to output dir --------------------
for f in output/spectra/bw_pi.dat output/spectra/bw_K.dat output/spectra/bw_p.dat; do
    [ -f "$f" ] && cp "$f" "$OUTPUT_DIR/spectra/"
done

echo ""
echo "=================================================="
echo " Done."
echo " Input :  $INPUT_DIR"
echo " Output:  $OUTPUT_DIR"
echo "   spectra/   phenix_raw.txt, phenix_parsed.csv, bw_*.dat"
echo "   figures/   phenix_spectra.png"
echo "   logs/      run.log"
echo "=================================================="