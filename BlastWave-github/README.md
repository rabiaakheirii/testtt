# BlastWave

## macros/

### `run_fit.C`

Main entry point: load a JSON config, run the selected fit mode, write parameters and PDFs.

```bash
root -l 'macros/run_fit.C("config/phenix_uu.json", "<mode>")'
```

Configs: `config/phenix_{uu,cuau,heau,pal,auau}.json`.

#### Fit modes

1. `load` — Load config and spectra only (no fit) |
2. `handNoFit` — Draw blast-wave curves from `hand_params` in JSON (no fitting) 
3. `handFit` — Individual fit per particle/centrality, `hand_params` use as init params
4. `global` — Hand fit, then global fit (common T, β for π+K+p per charge) 
5. `individualGlobal` — Global fit, then individual fits initialized from global results |
6. `full`  — Full chain: handFit → global → individualFromGlobal → global refine 

#### Draw modes (no refit)

1. `draw` — Spectra + individual curves from `BWparams.txt` 
2. `drawGlobal` — Spectra + global curves from `GlobalBWparams.txt` 
3. `drawParams` — Parameter plots: T(β), T(N<sub>part</sub>), β(N<sub>part</sub>) 

Optional `drawParams` filters (colon-separated):

```bash
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams:both:pos")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams:global")'
```

Sources: `individual` | `global` | `both`. Charge/species: `pos` | `neg` | `pi` | `k` | `p`. Add `nopredict` to hide the T(β) prediction curve.

### `draw_params_syst.C`

Overlay T(N<sub>part</sub>), β(N<sub>part</sub>), T(β) for several systems (reads existing fit txt files).

```bash
root -l 'macros/draw_params_syst.C("CuAu,UU", "global")'
```

### Output

- Parameters: `output/<system>/txtParams/BWparams.txt`, `GlobalBWparams.txt`
- PDFs: `output/<system>/`

---

## fit/

### `FitPipeline.h`

Orchestrates fit steps in order:

1. **HandNoFit** — set hand parameters, no Minuit  
2. **HandFit** — individual fit from hand init  
3. **Global** — simultaneous π+K+p fit (shared T, β)  
4. **IndividualFromGlobal** — individual fit starting from global results  
5. **GlobalRefine** — second global pass using updated individual params  

`RunFull()` runs steps 2–5. Results are written to txt via `ParamIO`.

### `IndividualFit.h`

Per-particle, per-centrality blast-wave fit (`IndividualFitter`).

- Init: hand params from JSON, or global params from `GlobalBWparams.txt`
- Limits: scaled around the init values, or absolute limits from config
- Strategy is a `FitStrategy { InitStrategy, LimitsStrategy }`

### `GlobalFit.h`

Global fit for one charge at a time (π, K, p with common T and β, separate normalizations). Uses Minuit2 / ROOT `Fit::Fitter`.

---

## config/

JSON configs per collision system: spectra path, centralities, N<sub>part</sub>, hand starting values, fit ranges/limits, draw options, output directory.

---

## core/

Shared infrastructure:

1. `FitConfig.h` — Load JSON into a typed config 
2. `SpectrumData.h` — Load p<sub>T</sub>/m<sub>T</sub> spectra 
3. `ParamIO.h` — Read/write individual and global parameter txt files 
4. `BlastWaveModel.h` — Blast-wave integrand / fit function 
5. `Constants.h` — Particle masses, names, array sizes 

---

## analysis/

Plotting helpers used by the macros:

1. `DrawSpectra.h` — Spectra PDFs with blast-wave curves (from fit files or handNoFit) 
2. `DrawParams.h` — Single-system T(β), T(N<sub>part</sub>), β(N<sub>part</sub>) 
3. `DrawParamsMulti.h` — Multi-system overlays of the same plots 

---

## Legacy

The following directories are **legacy** and are not part of the current workflow:

- **`PHENIX/`** — older standalone PHENIX macros and copies of spectra/params  
- **`mpdChargedHadrons/`** — earlier MPD / charged-hadron analysis macros  

Use `macros/`, `fit/`, `core/`, `config/`, and `analysis/` instead.
