# Theorem Inventory

## Validated Statements

### 1. The Core Identity
The MHV gravity amplitude (reconstructed from Hodges determinant or KLT) is proportional to the **Rooted Forest Polynomial** evaluated on the spinor edge variables $z_{ij}$.

\[ M_{MHV} = (-1)^{n-1} \frac{\langle 01 \rangle^8}{N_R \prod_{k \notin R} C_k^2} F_{n,R}(z_{ij}) \]

*   **Status:** Verified for $n=4, 5, 6$.
*   **Code:** `src/scripts/physics_pullback_n*.sage`.
*   **Relation to Literature:** Equivalent to NSVW tree formula + Matrix-Tree Theorem.

### 2. Forest Polytope Facets
The Newton polytope of $F_{n,R}$ (the "Forest Polytope") has boundaries that match physical factorization channels.
For $n=6$, the 26 facets classify into:
*   $x_{ij} \ge 0$ (14 facets): Contact terms / Simple boundaries.
*   $\sum_{e \in S} x_e \le k$ (subset sums): Correspond to soft/collinear limits where a subset of particles goes soft/collinear.
*   Global sum constraints.

*   **Status:** Verified for $n=6$.
*   **Code:** `src/posgeom/facets_to_subsets.py`.

### 3. Saddle Point Pushforward
The canonical form of the Forest Polytope $\Omega(P)$ can be computed as the pushforward of the dlog form on the positive orthant via the **Moment Map** $X = \nabla \log F(z)$.
\[ \Omega(X) = \sum_{z^* : X(z^*)=X} \frac{1}{\det J(z^*)} \]
This resolves the dimensional mismatch between the parameter space (108 dims for n=6) and the intrinsic polytope space (11 dims).

*   **Status:** Verified numerically for $n=4, 5$.
*   **Code:** `src/posgeom/saddle_pushforward.py`, `src/tests/test_saddle_pushforward.py`.

## Open Conjectures / Next Steps (Phase V)

1.  **Exact Jacobian Factor:** The map sweep for $n=6$ shows that the raw pushforward value is not simply equal to the amplitude (Ratio ~ $10^{-20}$ or $NaN$). There is likely a missing **Jacobian factor** or prefactor in the definition of the map itself to match the physical normalization.
2.  **Map Uniqueness:** Is the Moment Map the *unique* map that generates the correct geometry? Or is there a "Twisted" map?
3.  **Residues:** Explicitly matching the residues of the pushforward form to the factorization limits of gravity.
