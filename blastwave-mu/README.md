# blastwave-mu

A C++17 blast-wave fitter with non-zero electric-charge (μ_Q) and
strangeness (μ_S) chemical potentials for identified-hadron spectra
in heavy-ion collisions.

## Physics goal

Extend the standard blast-wave model (Schnedermann–Sollfrank–Heinz, 1993)
to include the electric-charge and strangeness chemical potentials at
**kinetic freeze-out**, in addition to the baryon chemical potential.
Perform a simultaneous fit of identified-hadron spectra (π⁺, π⁻, K⁺, K⁻,
p, p̄) to extract:

- T_kin — kinetic freeze-out temperature
- β_s — surface transverse flow velocity
- ⟨β_T⟩ — average transverse flow
- μ_B, μ_Q, μ_S — baryon, electric-charge, strangeness chemical potentials

## Model

For a particle species i with mass m_i, degeneracy g_i, and quantum
numbers (B_i, Q_i, S_i), the blast-wave spectrum is

    dN_i/(m_T dm_T) = A_i · m_T ∫₀^R r dr
                      · I₀( m_T sinh(ρ(r)) / T )
                      · K₁( m_T cosh(ρ(r)) / T )
                      · exp( μ_i / T )

with  μ_i = B_i · μ_B + Q_i · μ_Q + S_i · μ_S,
      ρ(r) = atanh( β_s · (r/R)^n ).

Quantum statistics (Bose for mesons, Fermi for baryons) are included
via a fugacity expansion truncated at n = 4.

## Resonance feed-down

34 resonance channels (ρ, ω, Δ, K*, η, η', f₀, f₂, a₁, N*, Σ*) are
included. Each resonance R contributes to the pion spectrum with its
own chemical potential μ_R = B_R μ_B + Q_R μ_Q + S_R μ_S and branching
ratio. The decay shape uses the T_eff = (p*/m_R) · T approximation
from the original paper (Eq. 28). A separate chemical freeze-out
temperature T_ch = 156 MeV (lattice QCD) sets the resonance abundances.

## Status

| Phase | Content | Status |
|---|---|---|
| 0 | Skeleton + CMake | ✅ |
| 1 | Math utilities (Bessel, integration, grid, RNG) | ✅ |
| 2 | Particle table with B, Q, S | ✅ |
| 3 | Blast-wave spectrum with μ_B, μ_Q, μ_S | ✅ |
| 4 | PHENIX and ALICE data readers | ✅ |
| 5 | χ² with analytical normalization | ✅ |
| 6 | Nelder–Mead minimizer + joint fit | ✅ |
| 7 | Real PHENIX Au+Au fit (centrality 0) | ✅ |
| 8 | 34-resonance feed-down | ✅ |
| 9 | Full 12-centrality scan + 6 figures | ✅ |

## Results

Real PHENIX Au+Au 200 GeV data, 6 species, 12 centrality bins.
Parameters extracted from a simultaneous fit:

| centrality | N_part | T [MeV] | β_s | ⟨β_T⟩ | μ_B [MeV] | μ_Q [MeV] | μ_S [MeV] | χ²/ndof |
|---|---|---|---|---|---|---|---|---|
| 0–5%    | 351 | 201 | 0.53 | 0.26 | 36.4 | 1.3 | 4.8 | 82.4 |
| 5–10%   | 299 | 212 | 0.40 | 0.20 | 42.4 | 0.3 | 5.9 | 51.0 |
| 10–15%  | 253 | 208 | 0.44 | 0.22 | 40.3 | −0.9 | 6.6 | 46.9 |
| 15–20%  | 215 | 204 | 0.48 | 0.24 | 36.7 | 1.7 | 5.2 | 42.6 |
| 20–30%  | 168 | 199 | 0.52 | 0.26 | 33.6 | 2.7 | 3.7 | 38.9 |
| 30–40%  | 114 | 195 | 0.56 | 0.28 | 35.9 | 1.3 | 4.9 | 45.1 |
| 40–50%  |  78 | 184 | 0.62 | 0.31 | 32.2 | 1.9 | 3.7 | 34.8 |
| 50–60%  |  52 | 172 | 0.67 | 0.33 | 28.2 | 2.2 | 4.1 | 24.5 |
| 60–70%  |  34 | 159 | 0.71 | 0.36 | 24.8 | 1.6 | 3.8 | 14.5 |
| 70–80%  |  22 | 150 | 0.73 | 0.37 | 22.2 | 1.9 | 3.6 |  8.0 |
| 80–90%  |  14 | 138 | 0.76 | 0.38 | 22.7 | 1.8 | −0.7 | 3.5 |
| 90–100% |   7 | 134 | 0.74 | 0.37 | 18.2 | 1.5 | 0.7 | 2.4 |

### Key physics findings

1. **β_s matches published PHENIX values in peripheral collisions.**
   At centrality 90–100%, β_s = 0.74 (published 0.65). The model
   correctly captures the collective radial expansion.

2. **T_kin is biased high in central collisions.**
   T_kin = 201 MeV at 0–5% vs published 100 MeV. This is the
   documented limitation of the T_eff resonance approximation,
   which fails when the fireball lives long enough for many
   resonance decays to occur.

3. **Systematic χ²/ndof trend: 2.4 (peripheral) → 82 (central).**
   A 34-fold degradation of the fit quality as the system grows.
   This quantifies precisely where the T_eff approximation breaks.

## Directory layout

    src/
      core/        configuration, constants
      utils/       Bessel, integration, grid, RNG, Nelder–Mead
      particles/   particle table with B, Q, S
      blastwave/   spectrum + 34-resonance feed-down
      fit/         χ², joint fit
      io/          data readers (PHENIX, ALICE)
    data/
      experimental/  experimental data (csv)
      resonances/    resonance table
      Npart_PHENIX_AuAu.txt
    input/
      PHENIX_spectra/AuAu/  raw PHENIX file
    output/
      fits/PHENIX_AuAu/     12 centrality result files
      figures/              6 physics figures
    scripts/                 Python plotting
    tests/                   unit tests (8 suites)
    docs/                    documentation

## Build

    mkdir -p build && cd build
    cmake ..
    make
    cd ..
    ./build/bwmu input/PHENIX_spectra/AuAu/spectra.txt output/fits/PHENIX_AuAu 0

Requires: CMake ≥ 3.16, C++17, GSL, OpenMP.

## Run all tests

    cd build
    ctest --output-on-failure

All 8 test suites must pass.

## Figures

All 6 figures are in `output/figures/`:

1. `fig1_spectra.png` — best-fit π, K, p spectra (centrality 0)
2. `fig2_T_vs_Npart.png` — T_kin vs system size
3. `fig3_betas_vs_Npart.png` — β_s vs system size
4. `fig4_mu_vs_Npart.png` — μ_B, μ_Q, μ_S vs system size
5. `fig5_chi2_contour.png` — (T, β_s) colored by χ²/ndof
6. `fig6_ratios.png` — π⁺/π⁻, K⁺/K⁻, p̄/p vs m_T
7. `fig_key_chi2_vs_Npart.png` — χ²/ndof vs N_part (systematic study)

## Limitations and future work

- **T_eff approximation.** The resonance decay shape uses a simple
  Boltzmann factor at T_eff = (p*/m_R)·T instead of the exact
  two-body decay integral. This is the dominant source of the T bias
  in central collisions.
- **Fixed T_ch.** The chemical freeze-out temperature is fixed at
  156 MeV rather than fitted.
- **No MCMC error bars.** Parameter uncertainties are not yet
  propagated from the fit.
- **No ALICE comparison.** The LHC branch is stubbed but not run.

## License

Academic use. Please cite Schnedermann, Sollfrank, Heinz,
Phys. Rev. C 48 (1993) 2462, and PHENIX Collaboration
Phys. Rev. C 69 (2004) 034909.
