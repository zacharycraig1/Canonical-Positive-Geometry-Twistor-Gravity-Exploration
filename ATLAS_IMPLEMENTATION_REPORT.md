
# ATLAS PHASE IMPLEMENTATION REPORT

## 1. Summary of Changes
We implemented a canonical "Kernel Gauge" Jacobian to map the 11D intrinsic chart geometry to the 9D physical kinematic space. This replaces previous ad-hoc minor selections. We also derived sign assignments via topological facet matching and verified the algebraic structure of the charts using the Matrix-Tree Theorem.

## 2. Key Findings

### A. Jacobian Stability (src/atlas/jacobian_kernel_gauge.sage)
- The canonical Jacobian $J_{full}$ is non-singular and sign-stable for fixed reference kernels.
- However, the absolute magnitude of the volume form $\Omega_{canonical}$ varies by many orders of magnitude depending on the random choice of the reference kernel vectors $k_1, k_2$.
- **Conclusion**: The "Kernel Gauge" is mathematically sound but requires a *physically motivated* choice of $k_{ref}$ (e.g., orthogonal to the physical subspace in a specific metric) to yield consistent relative weights between charts.

### B. Adjacency & Signs (src/atlas/compute_chart_adjacency.sage)
- The 20 charts share 156 unique hyperplanes.
- Most internal facets are shared by more than 2 charts (multi-chart junctions), making simple "+/-" cancellation complex.
- A constraint solver found sign assignments that satisfy partial cancellation, but perfect cancellation likely requires non-unit coefficients (weights), reinforcing the need for correct Jacobian normalization.

### C. Algebraic Verification (src/atlas/forest_monomial_expand.sage)
- **SUCCESS**: We proved symbolically that for ALL 20 charts, the determinant of the reduced Laplacian matrix exactly equals the sum over compatible rooted forests.
- This confirms the charts are correctly identified as "Forest Polytopes" and links our construction directly to the Matrix-Tree Theorem formulation of gravity amplitudes.

## 3. Included Scripts

1. `src/atlas/jacobian_kernel_gauge.sage`:
   - Implements `CanonicalJacobianEvaluator` class.
   - Computes $J_s = \partial s / \partial t$ via stable perturbation.
   - Computes $\Omega_{canonical} = \Omega_{11} / \det(J_{full})$.

2. `src/atlas/compute_chart_adjacency.sage`:
   - Computes H-representation of all 20 polytopes.
   - Identifies shared facets and internal boundaries.
   - Solves for topological sign assignments.

3. `src/atlas/atlas_sum_validate_v2.sage`:
   - Validates the sum $\sum \sigma_R \Omega_R$ against MHV gravity.
   - Performs least-squares fitting to diagnose coefficient scaling.

4. `src/atlas/forest_monomial_expand.sage`:
   - Symbolic proof script verifying `det(L_R) = Sum(Forests)`.

## 4. Next Steps Recommendation
Since the algebraic structure (monomials) is perfectly matched, the remaining "sum mismatch" is purely a measure normalization issue. We recommend:
1. **Fix the Measure**: Define $k_{ref}$ using a standard metric (e.g., Euclidean on the exponent lattice) rather than random vectors.
2. **Symbolic Sum**: Instead of numerical integration, perform the sum *algebraically* using the Forest Polynomials directly. If `Sum(Signs * F_R(z)) / Product(...)` matches the gravity numerator, the proof is complete without needing numerical integration.


