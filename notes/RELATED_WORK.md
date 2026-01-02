# Related Work and Prior Art

This document tracks literature relevant to the Phase U "Stringy Canonical Forms" and "Saddle Pushforward" approach for gravity amplitudes.

## Key References

### Stringy Canonical Forms & Positive Geometry
*   **Arkani-Hamed, He, Lam (2019)**: "Stringy Canonical Forms".
    *   Defines the stringy integral over the positive orthant $\int \prod \frac{dx}{x} x^{\alpha' X} p(x)^{-\alpha'}$.
    *   Establishes the link to polytope canonical forms in the $\alpha' \to 0$ limit.
    *   Discusses the "scattering equations" (saddle points) as the method to evaluate these integrals, involving the Jacobian determinant.
    *   **Relevance:** This is the mathematical framework we are implementing in Phase U. Our "Candidate A" map $X = \nabla \log p(z)$ is the moment map discussed here.

### MHV Gravity and Matrix-Tree Theorems
*   **Nguyen, Spradlin, Volovich, Wen (2009)**: "The Tree Formula for MHV Gravity Amplitudes".
    *   First explicit formula expressing MHV gravity as a sum over spanning trees/forests.
    *   **Relevance:** The "Forest Polynomial" $F_{n,R}(z)$ we use is exactly this object (or its generalized forest version).

*   **Hodges (2011)**: "New expressions for gravitational scattering amplitudes".
    *   Determinant formula for MHV gravity.
    *   **Relevance:** We use this as the "oracle" ground truth. The equivalence between the tree sum and the determinant is a consequence of the Matrix-Tree Theorem.

*   **Feng, He (2012)**: "KLT and Gravitational MHV Amplitudes".
    *   Explores the KLT relations and Matrix-Tree connections explicitly.
    *   **Relevance:** Connects the KLT kernel (used in our checks) to the Hodges determinant.

### Geometric Combinatorics
*   **Chaiken (1982)**: "A combinatorial proof of the all minors matrix tree theorem".
    *   **Relevance:** Provides the rigorous link between sums over rooted forests and minors of the Laplacian matrix (Hodges determinant).

*   **Postnikov (2009)**: "Permutohedra, Associahedra, and Beyond".
    *   Foundational work on polytopes related to graphs and matroids.
    *   **Relevance:** The Forest Polytope is a type of Generalized Permutohedron (or Polymatroid Base Polytope).

## Novelty Assessment (Phase U)
*   **Known:** The identity "MHV Gravity = Forest Polynomial(z)" is mathematically equivalent to the NSVW formula and Matrix-Tree theorem. It is **not new**.
*   **Plausibly New:** The construction of a **positive geometry** (via the Forest Polytope) and a **pushforward** map (via the Moment Map / Saddle Point) that recovers the amplitude.
    *   While AHBL discuss this for stringy integrals generally, applying it specifically to the *Forest Polytope of Gravity* with the *Spinor-Helicity Edge Variables* $z_{ij}$ as the parameter space is the specific contribution.
    *   The "Saddle Pushforward" mechanism resolving the dimension mismatch (108 parameters vs 16 intrinsic dimensions) via the Jacobian is the critical geometric insight.

## Theorem Inventory (Phase U Verified)
1.  **Combinatorial Identity:** $M_{MHV} \propto F_{n,R}(z_{ij})$ for n=4, 5, 6 (Verified).
2.  **Polytope Structure:** The Newton polytope of $F_{n,R}$ (Forest Polytope) has facets corresponding to physical soft/collinear limits (Verified for n=6).
3.  **Pushforward Mechanism:** The canonical form of the Forest Polytope is the $\alpha' \to 0$ limit of the stringy integral, and can be computed via sum over saddle points (Verified numerically for n=4, 5).
