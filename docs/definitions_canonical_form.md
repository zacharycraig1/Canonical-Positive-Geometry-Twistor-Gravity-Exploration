# Definitions: Positive Geometry & Canonical Forms

## 1. Positive Geometry
A **Positive Geometry** is a pair $(X, X_{\ge 0})$ consisting of a complex variety $X$ and a closed "positive part" $X_{\ge 0} \subset X(\mathbb{R})$, satisfying recursive boundary properties. 

In our context (MHV Gravity), the relevant geometry is the **3-Rooted Spanning Forest Polytope** $P_{n,R} \subset \mathbb{P}^N$.
- **X**: The Projective Space $\mathbb{P}^N$ (or a subspace/toric variety).
- **Positive Part**: The convex polytope $P_{n,R}$ itself.

## 2. Canonical Form
The **Canonical Form** $\Omega(X, X_{\ge 0})$ is a unique meromorphic differential form on $X$ characterized by:
1.  **Singularities**: It has logarithmic singularities (simple poles) exactly along the boundaries of $X_{\ge 0}$.
2.  **Residues**: The residue along a boundary component $C$ is the canonical form of the boundary geometry $\Omega(C, C_{\ge 0})$.
3.  **Normalization**: $\int_{X_{\ge 0}} \Omega = \pm 1$ (projective volume).

For a convex polytope $P \subset \mathbb{P}^d$ defined by vertices $Z_i$, the canonical form $\Omega(P)$ is a form on the **dual space** $(\mathbb{P}^d)^*$. If $W$ is a point in the dual space, the "canonical function" $\Omega(W)$ is the volume of the dual polytope.

### 2.1 Formula
For a polytope with vertices $Z_1, \dots, Z_k \in \mathbb{P}^d$, triangulated into simplices $S_a = (Z_{a_0}, \dots, Z_{a_d})$:
\[
\Omega(W) = \sum_{S_a} \text{sgn}(S_a) \frac{\det(Z_{a_0}, \dots, Z_{a_d})}{\prod_{i=0}^d (W \cdot Z_{a_i})}
\]

## 3. The Pushforward Definition
The main theorem of this work connects the physics (Gravity Amplitude) to the geometry via a **Pushforward**.

Let $X_{P}$ be the toric variety associated with the polytope $P_{n,R}$.
Let $\Omega(X_P)$ be the standard toric canonical form on $X_P$ (associated with the torus action).
Let $\phi: X_P \to (\mathbb{P}^N)^*$ be the "scattering map" (or algebraic moment map).

The Gravity Amplitude $M_n$ is the **pushforward** of the canonical form:
\[
M_n \, d\mu = \phi_* [\Omega(X_P)]
\]
Equivalently, the amplitude is the canonical function of the **dual** polytope $P^*$.

Our construction $F_{n,R}(z)$ is precisely the Newton Polynomial of the polytope $P$, and its reciprocal defines the canonical form on the dual.





