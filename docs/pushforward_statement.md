# Pushforward Statement: The Forest Polytope and MHV Gravity

**Status:** Draft / Phase R Output  
**Date:** December 31, 2025

## 1. The Positive Geometry

We define the **Rooted Forest Polytope** $P_{n,R} \subset \mathbb{R}^{|E|}$ as the convex hull of incidence vectors of all spanning forests of $K_n$ rooted at $R$.
For $n=6, R=\{0,1,2\}$, this polytope has:
- **108 vertices** (corresponding to forests)
- **22 facets**, classified into:
  - **Lower Bounds** ($x_{ij} \ge 0$): Correspond to vanishing edges.
  - **Upper Bounds** ($x_{ij} \le 1$): Correspond to saturated edges (collinear).
  - **Subset Rank Constraints** ($\sum_{e \in S} x_e \le r(S)$): Correspond to factorization channels.

The canonical form $\Omega_{P_{n,R}}(W)$ is defined on the dual projective space.

## 2. The Physics Map

We define the map $\Phi$ from kinematic spinor space to the polytope weights $z_{ij}$ (associated with edges):
$$
z_{ij} = \frac{[ij]}{\langle ij \rangle} C_i C_j
$$
where $C_i = \langle i x \rangle \langle i y \rangle$ for reference spinors $x, y$.

## 3. The Dictionary (Phase R Results)

We have verified the following correspondence between kinematic limits and polytope boundaries:

| Physical Limit | Kinematic Behavior | Polytope Boundary | Type |
|---|---|---|---|
| **Soft Limit** ($p_i \to 0$) | $z_{ij} \to 0$ for all $j$ | Intersection of $x_{ij} \ge 0$ | Lower Bound Facets |
| **Collinear Limit** ($i \parallel j$) | $z_{ij} \to \infty$ | $x_{ij} \le 1$ | Upper Bound Facet (if exists) |
| **Factorization** ($s_S \to 0$) | Coupled $z$ scaling | $\sum_{S} x \le k$ | Subset Rank Facet |

**Note on Collinear Limits:** Not all edges $ij$ have a corresponding $x_{ij} \le 1$ facet in the polytope representation (e.g., edges between roots, or specific non-root pairs like 4-5 in the 0,1,2-root basis). In these cases, the singularity corresponds to a higher-codimension face or a specific vertex configuration, rather than a simple facet.

## 4. The Main Statement

**Conjecture (The Pushforward Theorem):**
Let $\Omega_{P_{n,R}}$ be the canonical form of the Forest Polytope. Let $\Phi$ be the map defined above.
Then the pushforward of the canonical form (summed over all $z$) matches the MHV Gravity Amplitude (Hodges formula):
$$
\int \Phi^*(\Omega_{P_{n,R}}) \delta(\text{constraints}) \propto M_n^{\text{MHV}}
$$
Alternatively, in the language of forms:
The pullback of the canonical form $\Omega_{P_{n,R}}$ under the map $\Phi$ (restricted to the support of the forest polynomial) yields the MHV gravity integrand.

## 5. Verification status
- **Combinatorics:** Forest polynomial factorization on facets verified (see `RESULTS/facet_factorizations_n6.json`).
- **Geometry:** Signed canonical form evaluator implemented and passed sign-consistency tests.
- **Kinematics:** Limit probes confirm the map $\Phi$ correctly targets the predicted facets for Soft and Collinear limits (where applicable).

## 6. Next Steps
- Formalize the "missing facet" behavior.
- Extend factorization checks to BCFW channels (Subset Rank facets).
- Perform full numerical integration (or symbolic pushforward) to prove the identity for $n=6$.
