# QCD Collision Simulation

A pedagogical C++ implementation of the full heavy-ion collision
pipeline:

    percolation -> hydro (Bjorken) -> Cooper-Frye -> spectra -> CEP search

## Physics

- Phase 1-2: Thermal + resonance contribution to hadron mT spectra
  (Schnedermann-Sollfrank-Heinz, Phys. Rev. C 48 (1993) 2462).
- Phase 3:   Longitudinal flow with cutoff eta_max ~ 1.7.
- Phase 4:   Full blast-wave with transverse flow beta_s ~ 0.5.
- Phase 5:   String percolation initial conditions.
- Phase 6:   Relativistic hydrodynamics (1+1D Bjorken) with 3 EoS.
- Phase 7:   Cooper-Frye freeze-out.
- Phase 8:   Net-proton cumulants and critical point search.
- Phase 9:   Integration, summary, validation.

## Build

    mkdir build && cd build
    cmake ..
    make
    cd ..
    ./build/sim

Requires: CMake >= 3.16, GSL, OpenMP, C++17 compiler.

## Output

- output/spectra/     -- mT spectra, blast-wave fits, cumulants
- output/events/      -- percolation map, hydro profile
- output/figures/     -- all plots (scripts in scripts/)
- output/summary.txt  -- human-readable summary

## Git tags

phase0-done ... phase8-done, v1.0
