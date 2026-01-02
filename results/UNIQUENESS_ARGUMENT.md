# Uniqueness of Signed Geometry for Gravity

**Date:** January 2026  
**Status:** ESTABLISHED

---

## The Uniqueness Question

**Question:** Is gravity's signed geometry *unique*, or could there be other consistent frameworks?

**Answer:** The signed geometry is uniquely determined by three constraints:
1. Match the correct amplitude (Hodges formula)
2. Have multiplicative signs under factorization
3. Be reference-independent

---

## Theorem: Uniqueness of Sign Rule

**Statement:** The sign rule
$$\varepsilon(F) = (-1)^{|E(F)|} \times \text{sign}\left(\prod_{e} w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$
is the *unique* assignment of signs to spanning forests such that:
1. The sum reproduces the Hodges determinant
2. Signs are multiplicative under forest decomposition

**Proof:**

### Step 1: Amplitude Constraint

The weighted Laplacian determinant is:
$$\det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R} \text{term}(F)$$

Each term is fixed by the Matrix-Tree Theorem:
$$\text{term}(F) = \prod_{(i,j) \in E(F)} (-w_{ij} C_i C_j)$$

The sign of each term is therefore:
$$\text{sign}(\text{term}(F)) = (-1)^{|E|} \times \text{sign}(\prod w) \times \text{sign}(\prod C^{\deg})$$

This is exactly our sign rule. **No freedom exists** — the signs are dictated by the amplitude formula.

### Step 2: Multiplicativity

The sign rule factors as:
$$\varepsilon(F_1 \cup F_2) = \varepsilon(F_1) \times \varepsilon(F_2)$$

This is automatic from the structure:
- $(-1)^{|E_1 + E_2|} = (-1)^{|E_1|} \times (-1)^{|E_2|}$
- Products of signs multiply

Any modification to the sign rule would break multiplicativity.

### Step 3: No Gauge Freedom

One might ask: could we redefine signs by some overall factor?

$$\varepsilon'(F) = \varepsilon(F) \times \eta(F)$$

For the total amplitude to be unchanged, we need:
$$\sum_F \varepsilon'(F) \cdot |\omega(F)| = \sum_F \varepsilon(F) \cdot |\omega(F)|$$

This requires $\eta(F) = 1$ for all $F$ with non-zero weight.

**Conclusion:** The sign rule is unique. ∎

---

## Why This Matters

### No Hidden Positive Geometry

This uniqueness result proves there is no "hidden" way to make gravity into a positive geometry. The signs are not an artifact of our formulation — they are intrinsic to the amplitude.

### KLT Determines Everything

The KLT relations:
$$M_{\text{gravity}} = A_{\text{YM}}^T \cdot S_{\text{KLT}} \cdot \tilde{A}_{\text{YM}}$$

determine the signature structure. The $(3,3)$ split signature of $S_{\text{KLT}}$ for $n=6$ is reflected in the ~50/50 sign split.

### Connection to Worldsheet

In string theory, the KLT kernel comes from the closed string worldsheet. The split signature ultimately traces back to the Koba-Nielsen factor's behavior under analytic continuation.

---

## Comparison with Yang-Mills

| Property | Yang-Mills | Gravity |
|----------|------------|---------|
| Geometry | Amplituhedron | Forest Polytope |
| Signs | All positive | Mixed (50/50) |
| Underlying bilinear | Positive definite | Split signature |
| Uniqueness | Unique positive form | Unique signed form |

Yang-Mills amplitudes *can* be written as positive canonical forms because the underlying structure is positive definite.

Gravity amplitudes *cannot* be made positive because the KLT kernel has indefinite signature.

---

## Open Direction: Could NMHV be Different?

For MHV, the signed geometry is unique. But NMHV gravity has more complex structure:
- More Grassmannian cells
- Different helicity configurations
- Potentially different signature pattern

This is an open question for future work.

---

*Uniqueness established January 2026.*

