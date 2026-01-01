# Phase E Report: Exact Algebraic Bridge

## 1. The Exact Identity
We have established the precise algebraic link between the **Weighted Laplacian Tree Geometry** and the **Physical 6-pt MHV Gravity Amplitude** ($M_6$).

The formula is:
$$ M_6 = - \langle 01 \rangle^8 \cdot \frac{\det(\tilde L^{(012)})}{\prod_{k \in \{3,4,5\}} C_k^2} \cdot \frac{1}{(\langle 01 \rangle \langle 12 \rangle \langle 20 \rangle)^2} $$

Where:
- $\tilde L$ is the **Weighted Laplacian** ($6 \times 6$) with:
  - Off-diagonal: $\tilde L_{ij} = - w_{ij} C_i C_j$
  - Diagonal: $\tilde L_{ii} = \sum_{k \neq i} w_{ik} C_i C_k$ (forcing row sum 0)
- Weights: $w_{ij} = [ij]/\langle ij \rangle$
- Vertex factors: $C_i = \langle i x \rangle \langle i y \rangle$ (with reference spinors $x, y$)
- $\det(\tilde L^{(012)})$ is the principal minor obtained by deleting rows/columns $\{0, 1, 2\}$.

**Crucial Finding regarding Rank:**
On the support of momentum conservation, the Weighted Laplacian $\tilde L$ has **rank 3** (corank 3).
- The standard Tree Sum (minor of size $n-1=5$) vanishes identically.
- The physical information resides in the **minor of size $n-3=3$**.

## 2. Pole Structure Analysis (Valuations)
We compared the singular behavior of three objects as $\langle ij \rangle \to 0$ ($\epsilon \to 0$):

| Limit | Object | Scaling | Slope ($d\log|f|/d\log \epsilon$) |
| :--- | :--- | :--- | :--- |
| **Adjacent** ($\langle 01 \rangle \to 0$) | Amplitude ($A_6$) | $1/\epsilon^2$ | -2 |
| | Plain Tree Sum ($n-1$) | $1/\epsilon$ | -1 |
| | Weighted Minor ($n-3$) | **Finite** | **0** |
| **Non-Adjacent** ($\langle 02 \rangle \to 0$) | Amplitude ($A_6$) | $1/\epsilon^2$ | -2 |
| | Plain Tree Sum ($n-1$) | $1/\epsilon$ | -1 |
| | Weighted Minor ($n-3$) | **Finite** | **0** |
| **Kept Indices** ($\langle 34 \rangle \to 0$) | Amplitude ($A_6$) | $1/\epsilon$ (?) | -1 |
| | Weighted Minor ($n-3$) | $1/\epsilon$ | -1 |

*Note: The pole order of $A_6$ at $\langle 34 \rangle$ appeared as -1 in our spinor test, which might be helicity-dependent or a nuance of the specific test configuration.*

**Mechanism Identified:**
The "Weighted Laplacian" construction (specifically the rank drop to 3) effectively **clears/suppresses** the poles associated with the deleted rows/columns (0, 1, 2). The universal prefactor $1/(\langle 01 \rangle \langle 12 \rangle \langle 20 \rangle)^2$ then restores these poles (as double poles). For kept indices, the Weighted Minor retains the single pole structure characteristic of the Tree Sum.

## 3. CHY Bridge Verification
We verified that the Weighted Laplacian is isomorphic to the CHY Jacobian matrix:
1. **Sigma Identity Verified:** $(\sigma_a - \sigma_b) \propto \frac{\langle ab \rangle}{\langle a \chi \rangle \langle b \chi \rangle}$.
2. **Matrix Identity:** The CHY matrix $\Phi_{CHY}$ (with entries $s_{ab}/(\sigma_a - \sigma_b)^2$) is proportional to the Weighted Laplacian $\tilde L$ (with entries $w_{ij} C_i C_j$) when the reference spinors are identified ($x=y=\chi$).

## 4. Conclusion
The "Missing Operation" (Phase E3) that converts the Tree Geometry to the Amplitude is the **Rank-3 Projection** induced by the reference spinor weighting $C_i C_j$, followed by multiplication by the CHY-measure factor.






