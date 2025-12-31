# Master Results: The Positive Geometry of MHV Gravity

## 1. The Core Identity
The $n$-point MHV Gravity Amplitude is exactly the "Reference-Weighted Laplacian Minor", normalized by the CHY measure.

### Formula
$$ \mathcal{M}_n = (-1)^{n-1} \langle a b \rangle^8 \cdot \frac{\det(\tilde{L}^{(R)})}{\prod_{k \notin R} C_k^2} \cdot \frac{1}{|\mathcal{N}(R)|^2} $$

Where:
- $R$ is any set of 3 "root" indices (e.g. $\{0, 1, 2\}$).
- $a, b$ are the negative helicity legs.
- $\tilde{L}$ is the **Weighted Laplacian** with weights $w_{ij} = [ij]/\langle ij \rangle$ and vertex factors $C_i = \langle i x \rangle \langle i y \rangle$.
- $\mathcal{N}(R)$ is the normalization factor (product of angle brackets cyclic in $R$).

### Status
- **Verified for n=6:** Exact rational arithmetic (Phase E).
- **Verified for n=7:** Exact rational arithmetic (Phase F).
- **Gauge Invariance:** Verified independence of reference spinors $x, y$.

## 2. Geometric Interpretation: Rooted Forests
The determinant $\det(\tilde{L}^{(R)})$ expands into a sum over **spanning forests** of the complete graph $K_n$, consisting of exactly 3 trees, where each tree contains exactly one root from $R$.

$$ \det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R(K_n)} \prod_{(i,j) \in E(F)} (-w_{ij} C_i C_j) $$

- **Positivity:** On the positive kinematic region (moment curve), all $w_{ij} > 0$. The signs in the Laplacian expansion interact with the determinant sign to yield a uniform sign object (for appropriate choices of roots).
- **Polytope:** The Newton polytope of this object is the **3-Rooted Spanning Forest Polytope**, a generalized permutohedron of dimension $n-1$ (embedded in edge space). For $n=6, R=\{0,1,2\}$, it has 108 vertices.

## 3. Pole Mechanism
The Weighted Laplacian has **corank 3** on the support of momentum conservation.
- This rank deficiency "clears" the singularities associated with the deleted rows/columns (e.g., $\langle 01 \rangle \to 0$ leads to a finite minor).
- The universal prefactor restores these singularities as double poles $\langle 01 \rangle^{-2}$.
- Singularities associated with "kept" rows/columns (e.g. $\langle 34 \rangle \to 0$) manifest as simple poles in the minor, matching the amplitude scaling.

## 4. Connection to CHY
The Weighted Laplacian is isomorphic to the **CHY Jacobian Matrix**:
$$ \Phi_{CHY} \propto \tilde{L} $$
when the CHY auxiliary spinors match the Laplacian reference spinors.
- The identity $(\sigma_a - \sigma_b) \propto \frac{\langle ab \rangle}{\langle a \chi \rangle \langle b \chi \rangle}$ is the key link.

## 5. Next Steps
- **n=8+:** Numerical checks (determinant based).
- **Permutation Sums:** Exploring if the full permutation symmetric object simplifies further.
- **Loop Level:** Investigation of whether this Laplacian structure persists at 1-loop.



