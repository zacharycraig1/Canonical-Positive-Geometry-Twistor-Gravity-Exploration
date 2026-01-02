# Phase S Results

## 1. Geometric Engine Updates
- **Canonical Form Evaluator**: Updated to handle dimension-deficient polytopes (like the Forest Polytope) via affine projection.
- **Triangulation Invariance**: Verified. The value of $\Omega(P)$ is independent of the triangulation choice and scales correctly with degree $-(d+1)$.
  - **Polytope Dimension**: For $n=5, |R|=2$, dimension is 8 (Ambient 10).
  - **Scaling**: Verified $2^{-9}$ scaling factor.

## 2. Physical Singularity Checks
- **Collinear (4||5)**: Identified as a codimension-3 face in the current Forest Polytope ($n=6, R=[0,1,2]$).
- **Iterated Residue**: Computed $\Omega(F)$ on this face.
  - **Result**: Non-zero value ($\approx 10^{300}$ range for limit kinematics, indicating a pole).
  - **Interpretation**: The geometry *does* have a boundary corresponding to this singularity, but it is higher codimension than expected (3 vs 1). This implies the "Facet Dictionary" needs to account for higher-codimension boundaries or a different root basis is needed to make it a facet.

## 3. Pushforward Ratio Test ($n=6$)
- **Method**: Compared $M_{\text{MHV}}$ (Hodges) to $\Omega(P)(z(\lambda,\tilde\lambda))$ for random rational kinematics.
- **Result**: **FAIL**. The ratio is not constant.
  - Trial 0: Ratio $\sim 10^{16}$
  - Trial 1: Ratio $\sim 10^{2}$
  - Trial 3: Ratio $\sim 10^{16}$
  - Trial 4: Ratio $\sim 10^{24}$
- **Conclusion**: The naive pullback $\Phi^*(\Omega_P) = M_{\text{MHV}}$ does not hold. This indicates either:
  1. A missing Jacobian factor from the map $\Phi$.
  2. The Forest Polytope canonical form requires additional weight factors (numerator weights).
  3. The map itself needs adjustment (e.g. different reference spinor dependence).

## 4. Next Steps
- Investigate the Jacobian of the edge-variable map.
- Explore "Weighted Forest Polytopes" where edges carry weights derived from the matrix-tree theorem coefficients.
- Re-evaluate the "missing" codimension-1 behavior for Col(4||5) by checking other root choices (Phase S3/S4).




