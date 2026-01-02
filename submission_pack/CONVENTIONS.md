# Sign and Normalization Conventions

**Date:** January 2026  
**Purpose:** Document all sign/normalization conventions used in this project to ensure reproducibility and referee clarity.

---

## 1. Weighted Laplacian Definition

The Hodges weighted Laplacian is defined as:

$$\tilde{L}_{ij} = \begin{cases}
-w_{ij} \cdot C_i \cdot C_j & \text{if } i \neq j \\
\sum_{k \neq i} w_{ik} \cdot C_i \cdot C_k & \text{if } i = j
\end{cases}$$

where:
- $w_{ij} = [ij]/\langle ij \rangle$ is the kinematic edge weight
- $C_i = \langle i, x \rangle \langle i, y \rangle$ is the reference spinor factor
- $(x, y)$ are arbitrary reference spinors (the full amplitude is independent of this choice)

---

## 2. Matrix-Tree Theorem Conventions

### Standard MTT Convention

Define **positive edge weights**: $a_{ij} := -\tilde{L}_{ij} = w_{ij} C_i C_j$

Then by the All-Minors Matrix-Tree Theorem (Chaiken 1982):
$$\det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R} \prod_{e \in E(F)} a_e$$

### Signed-Edge Convention (Used in Code)

Our code computes forest terms using the actual Laplacian entries:
$$b_{ij} := \tilde{L}_{ij} = -a_{ij} = -w_{ij} C_i C_j$$

The relationship is:
$$\sum_{F} \prod_{e} b_e = (-1)^{|E|} \cdot \det(\tilde{L}^{(R)})$$

where $|E| = n - |R|$ is the number of edges per forest.

**Key files using this convention:**
- `src/signed_geometry/canonical_form.sage` (line ~125)
- `src/signed_geometry/verify_chy_sign_derivation.sage`

---

## 3. Sign Rule Formula

The sign of each forest term (using signed-edge convention) is:

$$\varepsilon(F) = (-1)^{|E(F)|} \times \text{sign}\left(\prod_{e} w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg_F(v)}\right)$$

| Factor | Origin | Dependence |
|--------|--------|------------|
| $(-1)^{\|E\|}$ | Convention choice (signed vs unsigned edges) | Combinatorial (constant for fixed n, R) |
| $\text{sign}(\prod w_e)$ | Kinematic ratios $[ij]/\langle ij\rangle$ | Kinematic (varies with momenta) |
| $\text{sign}(\prod C_v^{\deg})$ | Reference spinor factors | Gauge (cancels in full amplitude) |

---

## 4. Amplitude Formula

The full 6-point MHV gravity amplitude is:

$$\mathcal{M}_6 = (-1)^{n-1} \cdot \langle ab \rangle^8 \cdot \frac{\det(\tilde{L}^{(R)})}{\mathcal{N}_R \cdot \prod_{k \notin R} C_k^2}$$

where:
- $a, b$ are the two negative-helicity particles (typically 0, 1)
- $R = \{r_1, r_2, r_3\}$ is the root set (typically $\{0, 1, 2\}$)
- $\mathcal{N}_R = (\langle r_1 r_2 \rangle \langle r_2 r_3 \rangle \langle r_3 r_1 \rangle)^2$

The overall sign $(-1)^{n-1} = (-1)^5 = -1$ for $n=6$.

---

## 5. Spinor-Helicity Conventions

### Angle Bracket
$$\langle ij \rangle = \lambda_i^1 \lambda_j^2 - \lambda_i^2 \lambda_j^1 = \det(\lambda_i, \lambda_j)$$

### Square Bracket
$$[ij] = \tilde{\lambda}_i^1 \tilde{\lambda}_j^2 - \tilde{\lambda}_i^2 \tilde{\lambda}_j^1 = \det(\tilde{\lambda}_i, \tilde{\lambda}_j)$$

### Momentum Conservation
$$\sum_{i=1}^n \lambda_i^a \tilde{\lambda}_i^{\dot{a}} = 0$$

---

## 6. KLT Kernel Convention

The KLT relation is:
$$M_{\text{gravity}} = \sum_{\alpha, \beta} A_{\text{YM}}(\alpha) \cdot S[\alpha|\beta] \cdot \tilde{A}_{\text{YM}}(\beta)$$

The KLT kernel $S[\alpha|\beta]$ is computed using momentum kernel conventions from Bern-Carrasco-Johansson. For $n=6$, the kernel matrix has signature $(3, 3)$.

---

## 7. Code-to-Paper Mapping

| Paper Symbol | Code Variable | File |
|--------------|---------------|------|
| $\tilde{L}_{ij}$ | `L[i,j]` | `verify_chy_sign_derivation.sage` |
| $w_{ij}$ | `w[(i,j)]` | `verify_chy_sign_derivation.sage` |
| $C_i$ | `C[i]` | `verify_chy_sign_derivation.sage` |
| $\varepsilon(F)$ | `predicted_sign` | `compute_sign_from_factorization()` |
| Forest term | `compute_forest_term_direct()` | Uses signed-edge convention |

---

## 8. Conversion Between Conventions

To convert from signed-edge (code) to standard MTT (paper):

```python
# Code computes: sum_F prod_e (-w * C_i * C_j)
# Paper uses:    det(L) = sum_F prod_e (w * C_i * C_j)
# Relation:      code_sum = (-1)^|E| * det(L)

n, R = 6, 3
num_edges = n - len(R)  # = 3
det_L = (-1)**num_edges * code_forest_sum
```

---

*Conventions document created January 2026.*

