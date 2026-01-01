# Phase V Report: Exact Intrinsic Pushforward & Residues

**Date:** 2025-12-31
**Status:** Completed (Negative Result)

## Executive Summary
Phase V implemented the rigorous "Exact Intrinsic Basis" strategy to resolve normalization issues. While the geometric infrastructure (lattices, exact facets, moment map) is now numerically perfect, the **physics verification failed**.
1. The normalization ratio $M_{MHV} / \Omega$ is **not stable** across the moduli space (varies by $10^4$).
2. The **residue scaling** does not match: as $s_{01} \to 0$, Gravity scales as $1/\epsilon^2$, while Geometry $\Omega$ is constant ($\sim \epsilon^0$).

**Conclusion:** The current Forest Polytope + Map $z_{ij}$ combination does **not** capture the physical singularities of MHV gravity. The physical poles ($s_{ij} \to 0$) do not map to the geometric boundaries ($z_{e} \to 0$).

## 1. Geometric Achievements (Successes)
- **Facet Audit:** Resolved discrepancy. $n=6$ Forest Polytope has exactly **26 facets**. (`src/posgeom/facet_audit_n6.sage`)
- **Exact Lattice:** Implemented HNF-based integer basis, eliminating SVD errors. Covolume for $n=6$ is exactly $2\sqrt{3}$. (`src/posgeom/intrinsic_lattice.py`)
- **Moment Map:** Efficiently implemented using trace of Laplacian inverse. Verified $\det H = 0$ in ambient space but non-zero in projected space.

## 2. Physics Verification (Failures)
### A. Map Sweep Normalization
Run of `src/scripts/map_sweep_n6_exact.py` showed unstable ratios even with exact lattice:
- Sample 0: Ratio $\sim -1.8 \times 10^4$
- Sample 15: Ratio $\sim -1.0 \times 10^8$
- Normalized ratio (by $\langle xy \rangle^8$) also unstable.

### B. Residue Matching
Run of `src/posgeom/residue_match_n6.py` tested $s_{01} \to 0$ limit ($\epsilon \to 0$):
- **Gravity $M_{MHV}$**: Scales as $\epsilon^{-2.0}$ (Correct, pole).
- **Geometry $\Omega$**: Scales as $\epsilon^{-0.0008}$ (Constant, no pole).

This proves that the collinear limit $0 \to 1$ does not correspond to a boundary of the current Forest Polytope under the current map.

## 3. Next Steps
The failure is **structural**, not numerical.
- **Hypothesis:** The map $z_{ij}$ is insufficient or the polytope is "too small" (doesn't have the right boundaries).
- **Direction:** Re-evaluate the **physics map**. If $s_{ij} \to 0$ implies $z_{ij} \to 1$, we need a map where $z_{ij} \to 0$ or $\infty$.
- **Alternative:** The "Amplituhedron" for gravity might be the **Inverse Soft Factor** polytope or related to the **Hessian of the logarithm of the Amplitude** directly, rather than this specific Forest Polytope.



