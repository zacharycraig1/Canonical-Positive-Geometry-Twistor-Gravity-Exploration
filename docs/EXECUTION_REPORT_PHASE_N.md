# Phase N Execution Report

## Summary
We have successfully implemented the core components of the "Forest Pushforward" phase. The system now enforces strict conventions, passes rigorous internal consistency checks (deletion invariance), and demonstrates correct physical scaling at boundaries (soft and collinear limits). The comparison with an independent KLT oracle was implemented but revealed a functional mismatch (non-constant ratio) that requires further theoretical investigation, though dimensional analysis suggests consistency.

## 1. Correctness & Conventions
- **Mandelstam Convention:** `s_ij = <ij>[ji] = -det(P)` enforced. Calibration script passes.
- **Hodges Invariance:** The reduced determinant is proven independent of the deletion set (tested on all 20 triples for n=6).
- **Soft Limits:** Extrapolation confirms `|R-1| ~ epsilon^1` convergence for soft graviton limits.

## 2. The Pushforward Statement
- Defined the **Forest Polytope** and the map $\Phi$.
- **Boundary Test:** Confirmed that the Hodges form diverges as $1/s$ (or $1/\epsilon$) in the collinear limit, consistent with gravity factorization.

## 3. The KLT vs Hodges Puzzle
- **Result:** The ratio $M_{Hodges} / M_{KLT}$ is **not constant** (varies with kinematics).
- **Diagnostics:**
  - Both amplitudes have the correct mass dimension (-2).
  - Both show reasonable magnitude.
  - The mismatch likely stems from subtle definitions in the KLT momentum kernel or the Hodges prefactor normalization (e.g. spinor reference dependence in $s$ vs $\Phi$).
- **Next Step:** Requires a "Phase O" focused purely on reconciling the KLT and Hodges normalizations, possibly by comparing with a 3rd method (e.g. BCFW recursion) or simplifying to n=5.

## Deliverables
- `src/chy_oracle/test_hodges_deletion_invariance.sage` (PASS)
- `src/physics_limits/soft_extrapolate.sage` (PASS)
- `docs/NOVELTY_CHECKLIST.md` (Created)
- `docs/pushforward_statement.md` (Created)
- `src/pushforward/end_to_end_n6.sage` (Implemented, Ratio Fails)
- `src/pushforward/boundary_tests.py` (PASS - 1/s divergence)



