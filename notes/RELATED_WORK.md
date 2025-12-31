# Related Work & Positioning: Positive Geometry of Gravity

## Overview
This document clarifies the relationship between our "Positive Geometry" construction (Forest Polytope Canonical Form) and the existing literature on Gravity Amplitudes. The goal is to explicitly acknowledge known results to avoid "rediscovering" them, while highlighting our novel contributions regarding the geometric interpretation.

## 1. Known Results (Prior Art)

### The MHV Gravity Tree Formula
It is well-established that the MHV gravity amplitude can be expressed as a determinant of a "Hodges matrix" (or similar Laplacian-like matrices). 
- **Hodges (2012):** Provided the first determinant formula for MHV gravity amplitudes.
- **Nguyen et al. (2009):** Discussed tree formulas for gravity.

### Matrix-Tree Theorems
The mathematical identity linking determinants of Laplacian matrices to sums over spanning trees (or forests) is classical.
- **Matrix-Tree Theorem:** Relates the determinant of a reduced Laplacian to the sum over spanning trees.
- **Chaiken's All-Minors Matrix-Tree Theorem (1982):** Generalizes this to minors of any size, relating them to sums over rooted spanning forests.
- **Application to Physics:** It is known (though perhaps not always emphasized in geometric terms) that expanding the Hodges determinant yields a sum over trees/forests.

### Positive Geometries
- **Arkani-Hamed, Bai, Lam (2017):** Defined "Positive Geometries" and "Canonical Forms". They showed that the Amplituhedron and Associahedron are positive geometries for YM and bi-adjoint scalar theories, respectively.
- **Newton Polytopes & Stringy Canonical Forms:** The idea that a polynomial $F(z)$ defines a "stringy canonical form" whose singular locus is $F=0$ is part of the GKZ hypergeometric system framework.

## 2. Our Novelty Claim

We are **not** claiming the discovery that "Gravity MHV = Sum over Forests". This is mathematically equivalent to the Hodges formula via Chaiken's theorem.

**Our contributions are:**
1.  **Geometric Unification:** We identify the explicit **Positive Geometry** (the 3-rooted Spanning Forest Polytope $P_{n,R}$) whose *canonical form* generates the amplitude. This places Gravity on the same footing as the Amplituhedron (YM) and Associahedron (Scalars).
2.  **Pushforward Statement:** We formulate the amplitude as the **pushforward** of the canonical form on the toric variety $X_{P}$ associated with this polytope.
3.  **Factorization via Geometry:** We show that the *facets* of this specific polytope naturally encode the physical factorization channels of gravity (as double poles in $z_{ij}$).
4.  **Canonical Variables:** We identify the precise mapping from spinor kinematics to the toric "edge variables" $z_{ij}$ that makes this geometry manifest.

## 3. "Don't Rediscover" Checklist

| Topic | Likely Prior Art | What’s new in *your* work? | Status |
|---|---|---|---|
| Hodges determinant ↔ trees/forests | Hodges; Feng–He (matrix-tree) | Not new; we use this as a bridge to geometry. | Done |
| Positive geometry for gravity | Trnka/Paranjape “Gravituhedron” | We provide a concrete polytope construction for MHV | Done |
| Canonical forms via polytopes | ABL positive geometries; stringy | Identification of the *Forest Polytope* specifically | Done |
| Toric / Newton polytope pushforward | Stringy canonical forms | Application to the specific Gravity forest polynomial | Done |

## 4. Key References

*   **Hodges, A.** "A simple formula for gravitational MHV amplitudes." arXiv:1108.2227.
*   **Chaiken, S.** "A combinatorial proof of the all minors matrix tree theorem." SIAM J. Algebraic Discrete Methods, 1982.
*   **Arkani-Hamed, N., Bai, Y., Lam, T.** "Positive Geometries and Canonical Forms." arXiv:1703.04541.
*   **Bern, Z. et al.** "Gravity as the Square of Gauge Theory." (KLT Relations context).
*   **Feng, B., He, S.** "Graphs, determinants and gravity amplitudes."
