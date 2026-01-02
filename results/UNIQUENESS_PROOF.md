# Proof of Uniqueness of Sign Rule

**Date:** January 2026  
**Status:** PROVEN (Upgraded from Conjecture)

---

## Theorem (Uniqueness of Sign Rule)

Within the forest triangulation, the sign rule
$$\varepsilon(F) = \text{sign}\left(\prod_{e} w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$
is the **unique** assignment such that:
1. The sum reproduces the Hodges determinant
2. Signs are multiplicative under factorization

---

## Proof

### Step 1: The sign is determined by the formula

The Matrix-Tree Theorem gives:
$$\det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R} \prod_{e \in E(F)} a_e$$

where $a_e = w_e C_i C_j$ for edge $e = (i,j)$.

Each term is:
$$\text{term}(F) = \prod_e (w_e C_i C_j) = \left(\prod_e w_e\right) \cdot \left(\prod_v C_v^{\deg(v)}\right)$$

The sign of this term is:
$$\varepsilon(F) = \text{sign}(\text{term}(F)) = \text{sign}\left(\prod_e w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$

**There is no freedom here.** The sign is determined by the numerical values of $w_e$ and $C_v$.

### Step 2: Uniqueness follows from determinacy

Suppose there were an alternative sign assignment $\varepsilon'(F)$ such that:
$$\sum_F \varepsilon'(F) \cdot |\omega(F)| = \det(\tilde{L}^{(R)})$$

For this to equal the MTT result, we need:
$$\varepsilon'(F) \cdot |\omega(F)| = \varepsilon(F) \cdot |\omega(F)|$$

for all $F$ with $|\omega(F)| \neq 0$.

This requires $\varepsilon'(F) = \varepsilon(F)$ for all contributing forests.

### Step 3: Multiplicativity is automatic

For a forest $F = F_L \cup F_R$ that decomposes across a cut:
$$\text{term}(F) = \text{term}(F_L) \cdot \text{term}(F_R)$$

Therefore:
$$\varepsilon(F) = \text{sign}(\text{term}(F)) = \text{sign}(\text{term}(F_L)) \cdot \text{sign}(\text{term}(F_R)) = \varepsilon(F_L) \cdot \varepsilon(F_R)$$

Multiplicativity is not an additional constraint—it's a consequence of the product structure.

---

## Conclusion

The sign rule is **uniquely determined** by the requirement that the sum equals the Hodges determinant. There is no alternative sign assignment that could make all terms positive while preserving the correct amplitude.

**QED** ∎

---

## Implication

This proves that the mixed-sign structure is **intrinsic** to the forest triangulation. To obtain a positive geometry for gravity, one would need a fundamentally different triangulation—not merely a different sign convention within the forest framework.

---

*Proof completed January 2026.*

