# Phase B Update: Denominator Discovery Harness Implemented

## Status
- **Harness Refactored**: `src/scripts/discover_denominator.py` has been completely rewritten to follow the 3-stage robust workflow recommended in `chy_line_reconstruction_harness.md`.
- **Key Improvements**:
  1. **Stage A (Point Collection)**: Deterministic, single-threaded collection of exactly `target_pts` valid points. Swallow exceptions but log them. No silent partial failures.
  2. **Stage B (Rational Reconstruction)**: Adaptive degree search (checking total degree $D=0 \dots 100$) with holdout validation. Fixes $q_0=1$ to avoid ambiguity.
  3. **Stage C (Pole Analysis)**: **Removed QQbar root finding**. Now uses exact polynomial division (`GCD`) by precomputed quadratic $\langle ij \rangle(t)$ polynomials. This is numerically exact and much faster.

## Execution Issue
- Attempted to run the new harness, but **Docker is currently unreachable** (`npipe:////./pipe/dockerDesktopLinuxEngine` not found).
- The code is ready and optimized. It will require a working Sage/Docker environment to execute.

## Next Steps (Once Docker is up)
1. Run `src/scripts/discover_denominator.py`.
2. It will:
   - Collect 150 points.
   - Reconstruct $M_{CHY}(t) = P(t)/Q(t)$.
   - Factor $Q(t)$ into powers of $\langle ij \rangle(t)$.
   - Report the exact physical pole powers (e.g., $\langle ij \rangle^1$, $\langle ij \rangle^2$, etc.).
3. If successful, this uniquely identifies the clearing denominator $D(Z)$, allowing us to define the polynomial numerator $N(Z)$.



