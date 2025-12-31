# Theorem: The Positive Geometry of MHV Gravity Amplitudes

## 1. Geometric Object: The Weighted Laplacian
We define a **Weighted Laplacian** $\tilde{L} \in \mathbb{C}^{n \times n}$ associated with $n$ massless particles.

### Definitions
1.  **Kinematic Weights** ($w_{ij}$):
    For any pair of particles $i, j$, we define the sign-definite weight:
    $$ w_{ij} = \frac{[ij]}{\langle ij \rangle} $$
    On the "positive region" (bounded by the moment curve), $w_{ij} > 0$.

2.  **Reference Gauge** ($C_i$):
    Let $x, y \in \mathbb{CP}^1$ be arbitrary reference spinors. We define vertex weights:
    $$ C_i = \langle i x \rangle \langle i y \rangle $$

3.  **Weighted Laplacian** ($\tilde{L}$):
    $$ \tilde{L}_{ij} = \begin{cases} -w_{ij} C_i C_j & i \neq j \\ \sum_{k \neq i} w_{ik} C_i C_k & i = j \end{cases} $$
    By construction, $\sum_j \tilde{L}_{ij} = 0$ (row sums vanish).

## 2. The Master Identity
The $n$-point MHV Gravity Amplitude $\mathcal{M}_n$ is given by a specific **principal minor** of the Weighted Laplacian, corrected by the reference factors.

**Theorem:**
For any choice of reference spinors $x, y$ (generic), and any choice of 3 "deleted" indices $R=\{r, s, t\}$ (roots), the amplitude is:

$$ \mathcal{M}_n = (-1)^{n-1} \left( \prod_{a=1}^n \langle a, a+1 \rangle^{-2} \right) \times \frac{\det(\tilde{L}^{(R)})}{\prod_{k \notin R} C_k^2} \times \mathcal{J}_{\text{norm}}(R) $$

For the standard choice $R=\{0, 1, 2\}$ at $n=6$, with the standard Parke-Taylor-like prefactor normalization, this simplifies to our verified **Master Identity**:

$$ \mathcal{M}_6 = - \langle 01 \rangle^8 \cdot \frac{\det(\tilde L^{(012)})}{\prod_{k=3,4,5} C_k^2} \cdot \frac{1}{(\langle 01 \rangle \langle 12 \rangle \langle 20 \rangle)^2} $$

## 3. Corank 3 and Projections
Unlike the standard Laplacian (corank 1), the Weighted Laplacian $\tilde{L}$ has **corank 3** on the support of momentum conservation.
- This means all minors of size $>(n-3)$ vanish.
- The physical information is supported on the principal minors of size $(n-3) \times (n-3)$.
- This rank drop is the geometric mechanism that "clears" the poles associated with the deleted rows/columns.

## 4. Combinatorial Expansion (Forest Sum)
By the **All-Minors Matrix-Tree Theorem**, the determinant $\det(\tilde{L}^{(R)})$ expands into a sum over **spanning forests**:

$$ \det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R(K_n)} \prod_{(i,j) \in E(F)} (-w_{ij} C_i C_j) $$

where $\mathcal{F}_R(K_n)$ is the set of all spanning forests of the complete graph $K_n$ consisting of exactly 3 trees, where each tree contains exactly one root from $R=\{r, s, t\}$.

Substituting this back, the $C_i$ factors cancel in a specific way (to be proven gauge invariant), leaving a sum over positive-weight forests.

## 5. Gauge Invariance
The expression is invariant under the "gauge transformations" of the reference spinors:
- **Scaling:** $x \to \alpha x$, $y \to \beta y$.
    - $\tilde{L} \to (\alpha\beta)^2 \tilde{L}$
    - $\det(\tilde{L}^{(n-3)}) \to (\alpha\beta)^{2(n-3)} \det(\tilde{L}^{(n-3)})$
    - $\prod_{k \notin R} C_k^2 \to (\alpha\beta)^{2(n-3)} \prod C_k^2$
    - The ratio is invariant.
- **Shift:** $x \to x + \delta \lambda_i$ (checked numerically to be invariant).

This confirms the amplitude is a well-defined function of the kinematics $\lambda, \tilde\lambda$ alone, despite the intermediate introduction of $x, y$.

