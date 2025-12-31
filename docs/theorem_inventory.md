# Theorem Inventory: MHV Gravity & Positive Geometry

This document tracks the precise mathematical statements verified in this repository, distinguishing between **background knowledge** (known theorems) and **new claims** (novel geometric constructions).

## 1. Background Theorems (The "Ansatz")

These are established results we treat as infrastructure. Our contribution is **not** rediscovering these, but providing a specific geometric realization for them.

### T1. The Matrix-Tree Theorem for All Minors
**Source:** Chaiken (1982), "A combinatorial proof of the all minors matrix tree theorem".
**Statement:**
For a graph $G$ with Laplacian matrix $L$, the determinant of a submatrix $L^{(R)}_{(\bar{R})}$ (removing rows $R$ and columns $\bar{R}$) enumerates forests with specific root conditions.
**In our context:**
The "Weighted Laplacian" minor $\det(\tilde{L}^{(R)})$ exactly enumerates spanning forests where each component contains exactly one root from the set $R$.
\[
\det(\tilde{L}^{(R)}) = \sum_{F \in \text{Forests}(n, R)} \prod_{(i,j) \in E(F)} z_{ij}
\]
where $z_{ij}$ are the edge weights.

### T2. MHV Gravity as a Determinant
**Source:** Hodges (2011), "New expressions for gravitational scattering amplitudes".
**Statement:**
The $n$-point MHV gravity amplitude is given by:
\[
M_n^{\text{MHV}} = (-1)^{n-1} \sigma_{n,R} \det(\Phi^{(R)})
\]
where $\Phi$ is a Hodges matrix (weighted Laplacian).

### T3. Canonical Forms of Polytopes
**Source:** Arkani-Hamed, Bai, Lam (2017), "Positive Geometries and Canonical Forms".
**Statement:**
For a convex polytope $P$, there exists a unique canonical form $\Omega(P)$ on the dual projective space, which has logarithmic singularities on the boundaries dual to the vertices/facets of $P$.

---

## 2. New Geometric Claims (The "Result")

These are the specific statements we aim to prove and publish.

### C1. The "Forest Polytope" Identity
**Status:** Verified for $n=4, 5$.
**Statement:**
The Newton polytope of the MHV gravity numerator (in the specific edge variables $z_{ij}$) is the **3-rooted spanning forest polytope** of the complete graph $K_n$.
\[
P_{\text{grav}} = \text{Newton}\left( \det(\tilde{L}^{(R)}) \right) = \text{Conv}\left( \{ e_F \mid F \in \text{Forests}(n, R) \} \right)
\]

### C2. The Toric Pushforward (Target of Phase H)
**Status:** In Progress.
**Statement:**
There exists a projective toric variety $X_{P} \subset \mathbb{P}^{V-1}$ and a monomial map $\phi: (\mathbb{C}^*)^d \to X_P$ such that the pushforward of the standard toric form $\Omega_T = d\log t_1 \wedge \dots \wedge d\log t_d$ is the canonical form of the polytope $P$ (which equals the Amplitude).

**Precision required:**
We must define the map $\phi$ explicitly in terms of the lattice basis of the Forest Polytope.

### C3. Factorization via Facet Geometry
**Status:** To be verified.
**Statement:**
The facets of the Forest Polytope $P_{\text{grav}}$ are in one-to-one correspondence with the physical factorization channels of the amplitude (compatible with the poles of $\det(\tilde{L})$).



