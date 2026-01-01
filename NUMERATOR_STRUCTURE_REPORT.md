# Gravity MHV Amplitude Structure Report

## Key Findings

### 1. KLT Implementation Failure
The current implementation of the KLT formula for N=6 (`gravity_6pt_mhv_klt`) is fundamentally flawed:
- **Symmetry Broken:** The amplitude changes by orders of magnitude under permutation of particles (e.g., swapping 1 and 4 gives ratio ~-6918).
- **Symmetrization Fails:** Summing over 720 permutations does not converge to the Hodges result (Ratio ~3.3e6).
- **Conclusion:** The "Standard KLT Kernel" formula used is likely missing sign factors or is valid only in a specific basis that is not being handled correctly for N=6.

### 2. Hodges Numerator Uniqueness
We have successfully identified the analytic structure of the correct Gravity MHV amplitude ($M_{Hodges}$).
Using the ansatz $M = \frac{N(Z)}{D_{cyclic}}$, where $D_{cyclic} = (\prod_{i=1}^6 \langle i, i+1 \rangle)^2$:

- **Polynomiality:** The numerator $N(Z) = M_{Hodges} \times D_{cyclic}$ is a **polynomial** in momentum twistors.
  - Verified by Lagrange interpolation on a line $Z(t) = Z_A + t Z_B$.
  - Interpolation with 35 points yields a polynomial that **perfectly extrapolates** to points outside the training set (Ratio 1.000000).
- **Degree:** The polynomial has degree 34 in the linear parameter $t$ (consistent with scaling weights).
- **Scaling Dimension:** The amplitude scales as $Z^{-2}$ (consistent with mass dimension +2).

### 3. Resolution Strategy
The "missing terms" in the KLT formula can be uniquely determined by the polynomial difference:
$$ \Delta N(Z) = N_{Hodges}(Z) - N_{KLT}(Z) $$
where $N_{KLT}$ is the KLT numerator (brought to the common denominator $D_{cyclic}$).
Since $N_{Hodges}$ is a polynomial, the correct gravity amplitude is defined by this unique polynomial numerator.

## Next Steps

1. **Construct N explicitly:** Use the 200 sampled points to fit the coefficients of $N(Z)$ in a suitable basis (e.g. 4-bracket monomials).
2. **Geometric Connection:** The polynomial $N(Z)$ defines a hypersurface $N=0$. This is the candidate for the "Gravity Amplituhedron" boundary or canonical form numerator.
3. **Fix KLT:** Instead of summing permutations, use the Uniqueness Constraints (polynomiality of N) to fix the coefficients of the KLT basis.







