# Rigorous Proof of Uniqueness of Forest Coefficients

**Date:** January 2026  
**Status:** PROVEN (Polynomial-Monomial Approach)

---

## Theorem (Uniqueness of Monomial Coefficients)

In the polynomial ring $\mathbb{Z}[a_{ij}]$ with $a_{ij} = w_{ij} C_i C_j$ treated as independent formal variables, the forest expansion

$$\det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R} \prod_{e \in E(F)} a_e$$

has **unique** coefficients: each forest $F$ corresponds to a distinct monomial $\prod_{e \in E(F)} a_e$, and all coefficients are $+1$.

---

## Proof

### Step 1: Distinct Monomials

Different forests have different edge sets:
$$F \neq F' \implies E(F) \neq E(F')$$

Therefore their monomials are distinct:
$$\prod_{e \in E(F)} a_e \neq \prod_{e \in E(F')} a_e$$

in the polynomial ring $\mathbb{Z}[a_{ij}]$.

### Step 2: Unique Polynomial Expansion

A polynomial in independent variables has a **unique** expansion as a sum of monomials. Each monomial appears with a unique coefficient.

### Step 3: Apply Matrix-Tree Theorem

The All-Minors Matrix-Tree Theorem states:
$$\det(L^{(R)}) = \sum_{F \in \mathcal{F}_R} \prod_{e \in E(F)} a_e$$

where each forest contributes with coefficient $+1$.

### Step 4: Conclusion

Since the monomials are distinct and polynomials have unique expansions, the coefficient of each forest monomial is uniquely determined to be $+1$.

**QED** ∎

---

## Corollary: Signs on Real Slices

When we specialize to a **real kinematic slice** (rational spinors, split signature), each $a_e = w_e C_i C_j$ takes a real value.

The "sign of forest $F$" is:
$$\varepsilon(F) = \text{sign}\left(\prod_{e \in E(F)} a_e\right) = \text{sign}(\prod w) \times \text{sign}(\prod C^{\deg})$$

This sign is determined by the numerical values of $w_e$ and $C_v$ at the given kinematic point.

---

## What This Does and Does NOT Prove

### What It Proves:
- The forest expansion in the MTT basis has **unique, fixed coefficients**
- Any attempt to "re-sign" the terms changes the polynomial
- The mixed-sign structure is intrinsic to the forest triangulation

### What It Does NOT Prove:
- That no positive geometry exists for gravity
- That no other triangulation could yield uniform signs
- That the sign structure persists for complex kinematics

---

## Previous Incorrect Argument (Corrected)

The earlier version claimed: "If all terms had different signs, the sum would differ."

This is **not valid** without the polynomial structure. Sums can be equal even with different term-by-term values.

The correct argument uses **polynomial uniqueness**: distinct monomials have unique coefficients.

---

*Proof corrected January 2026.*
