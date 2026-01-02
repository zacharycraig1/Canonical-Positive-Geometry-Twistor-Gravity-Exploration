# CHY to Forest Polynomial Derivation

**Date:** January 2026  
**Status:** THEORETICAL DERIVATION

---

## Overview

This document traces how the CHY (Cachazo-He-Yuan) formula for MHV gravity amplitude produces the forest polynomial expansion, and identifies the origin of each sign factor in the discovered sign rule:

$$\varepsilon(F) = (-1)^{|E(F)|} \times \text{sign}\left(\prod_e w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$

---

## Step 1: The CHY Formula

The CHY formula for n-point gravity amplitude is:

$$\mathcal{M}_n^{\text{GR}} = \int d\mu_n \, (\text{Pf}'\Psi)^2$$

where:
- $d\mu_n$ = CHY measure localized on scattering equation solutions
- $\Psi$ = $(2n) \times (2n)$ antisymmetric matrix
- $\text{Pf}'\Psi$ = reduced Pfaffian (removing 2 rows/columns)

For **MHV amplitudes** (with particles $a, b$ having negative helicity), the formula simplifies:

$$\mathcal{M}_n^{\text{MHV}} = \langle ab \rangle^8 \times (\text{reduced determinant structure})$$

---

## Step 2: From Pfaffian to Determinant

For MHV gravity, the reduced Pfaffian squared relates to a determinant:

$$(\text{Pf}'\Psi)^2 = \det(\tilde{\Phi})$$

where $\tilde{\Phi}$ is the **Hodges matrix**:

$$\tilde{\Phi}_{ij} = \frac{[ij]}{\langle ij \rangle} = w_{ij} \quad (i \neq j)$$

The diagonal is determined by momentum conservation:

$$\tilde{\Phi}_{ii} = -\sum_{j \neq i} w_{ij}$$

This is a **weighted Laplacian** structure!

---

## Step 3: Weighted Laplacian and Matrix-Tree Theorem

The Hodges matrix has the structure of a graph Laplacian with weights $w_{ij} = [ij]/\langle ij \rangle$:

$$L_{ij} = \begin{cases} 
-w_{ij} & i \neq j \\
\sum_{k \neq i} w_{ik} & i = j
\end{cases}$$

The **generalized Matrix-Tree Theorem** states:

$$\det(L^{(R)}) = \sum_{F \in \mathcal{F}_R(K_n)} \prod_{(i,j) \in E(F)} w_{ij}$$

where:
- $L^{(R)}$ = minor with rows/columns in root set $R$ deleted
- $\mathcal{F}_R(K_n)$ = k-rooted spanning forests with root set $R$
- $|R| = k$ (for MHV, $k = 3$)

**This is the source of the forest sum!**

---

## Step 4: Including Vertex Factors (Reference Spinors)

The actual Hodges formula includes reference spinor-dependent factors:

$$\tilde{L}_{ij} = -w_{ij} \cdot C_i \cdot C_j \quad (i \neq j)$$

where $C_i = \langle i, x \rangle \langle i, y \rangle$ for reference spinors $x, y$.

The Matrix-Tree Theorem gives:

$$\det(\tilde{L}^{(R)}) = \sum_{F \in \mathcal{F}_R} \prod_{(i,j) \in E(F)} (-w_{ij} \cdot C_i \cdot C_j)$$

**Each forest term has the form:**

$$\text{term}(F) = \prod_{(i,j) \in E(F)} (-w_{ij} \cdot C_i \cdot C_j)$$

---

## Step 5: Factorizing the Forest Term

For a forest $F$ with edge set $E(F)$, let's factor the product:

$$\prod_{(i,j) \in E(F)} (-w_{ij} \cdot C_i \cdot C_j) = (-1)^{|E(F)|} \cdot \prod_{e \in E(F)} w_e \cdot \prod_{(i,j) \in E(F)} C_i \cdot C_j$$

Now, each vertex $v$ appears in the product $\prod C_i \cdot C_j$ exactly $\deg_F(v)$ times (once for each edge incident to $v$):

$$\prod_{(i,j) \in E(F)} C_i \cdot C_j = \prod_{v \in V} C_v^{\deg_F(v)}$$

Therefore:

$$\text{term}(F) = (-1)^{|E(F)|} \cdot \prod_{e \in E(F)} w_e \cdot \prod_{v \in V} C_v^{\deg_F(v)}$$

---

## Step 6: Origin of Each Sign Factor

### Factor 1: $(-1)^{|E(F)|}$

**Origin:** The Laplacian off-diagonal sign convention.

The weighted Laplacian has:
- $L_{ij} = -w_{ij} \cdot C_i \cdot C_j$ for $i \neq j$ (note the minus sign!)
- Matrix-Tree Theorem counts products of $(-L_{ij}) = w_{ij} \cdot C_i \cdot C_j$

But we use the "opposite" convention in the forest sum, so each edge contributes a factor of $(-1)$.

For $n = 6$ with $k = 3$ roots: $|E(F)| = n - k = 3$, giving $(-1)^3 = -1$.

---

### Factor 2: $\text{sign}(\prod_e w_e)$

**Origin:** Kinematic structure of spinor brackets.

The edge weight $w_{ij} = [ij]/\langle ij \rangle$ is a ratio of spinor brackets:
- $[ij]$ = square bracket (tilde spinors)
- $\langle ij \rangle$ = angle bracket (spinors)

The sign depends on:
1. The relative orientation of spinors $i$ and $j$
2. The kinematic configuration

**This sign is purely kinematic** - it doesn't depend on reference spinors or combinatorics.

For a given kinematic configuration, some edges have $w_{ij} > 0$, others have $w_{ij} < 0$.

---

### Factor 3: $\text{sign}(\prod_v C_v^{\deg_F(v)})$

**Origin:** Reference spinor choice in the Hodges formula.

The factor $C_v = \langle v, x \rangle \langle v, y \rangle$ depends on:
- The particle's spinor $\lambda_v$
- The reference spinors $x, y$ (gauge choice)

For a forest $F$, the product $\prod_v C_v^{\deg_F(v)}$ can be positive or negative depending on:
1. Which vertices have $C_v < 0$
2. The parities of vertex degrees $\deg_F(v)$

**Key insight:** While individual $C_v$ signs depend on reference choice, the **total amplitude is reference-independent** because the $C$ factors in the numerator cancel with $C$ factors in the denominator of the Hodges formula.

---

## Step 7: The Complete Sign Rule

Combining all factors:

$$\varepsilon(F) = (-1)^{|E(F)|} \times \text{sign}\left(\prod_{e \in E(F)} w_e\right) \times \text{sign}\left(\prod_{v \in V} C_v^{\deg_F(v)}\right)$$

And the weight magnitude:

$$\omega(F) = \left|\prod_{e \in E(F)} w_e \cdot \prod_v C_v^{\deg_F(v)}\right|$$

So each forest term is:

$$\text{term}(F) = \varepsilon(F) \cdot \omega(F)$$

---

## Step 8: Why ~50/50 Split?

The 54/54 split of positive and negative forest signs for $n = 6$ is **not coincidental**:

### Laplacian Sign Structure

The Laplacian has rank $n - 1$. For the reduced determinant (rank $n - 3$ for MHV), the cofactor expansion involves products of 3 off-diagonal elements, each with a $(-1)$ sign.

### KLT Kernel Connection

The KLT kernel $S_{\text{KLT}}$ relates gravity to YM amplitudes:
$$M_{\text{gravity}} = A_{\text{YM}}^T \cdot S_{\text{KLT}} \cdot \tilde{A}_{\text{YM}}$$

The kernel has **split signature (3,3)** for $n = 6$:
- 3 positive eigenvalues
- 3 negative eigenvalues

This split signature is the **geometric origin** of the 50/50 sign split!

### Twisted Cohomology

In the twisted cohomology framework:
$$S_{\text{KLT}} \propto (m_{\text{biadjoint}})^{-1}$$

The bi-adjoint matrix $m$ has the same split signature. The inverse preserves signature properties.

The split signature comes from the **twist potential** $W = \sum_{i<j} s_{ij} \log(z_i - z_j)$ having critical points with alternating Hessian signs.

---

## Step 9: Summary

| Sign Factor | Origin | Dependence |
|-------------|--------|------------|
| $(-1)^{\|E\|}$ | Laplacian off-diagonal convention | Combinatorial (fixed for given $n, k$) |
| $\text{sign}(\prod w)$ | Kinematic edge orientations | Kinematic only |
| $\text{sign}(\prod C^{\deg})$ | Reference spinor gauge choice | Reference + forest topology |

### The Full Picture

```
CHY Formula
    │
    ▼
Pfaffian² = Determinant
    │
    ▼
Hodges Matrix = Weighted Laplacian
    │
    ▼
Matrix-Tree Theorem → Forest Sum
    │
    ▼
Each term = (-w_ij C_i C_j) products
    │
    ▼
Sign Rule: ε(F) = (-1)^|E| × sign(∏w) × sign(∏C^deg)
```

---

## References

1. **CHY Formula:** Cachazo, He, Yuan - "Scattering of Massless Particles"
2. **Matrix-Tree Theorem:** Chaiken - "A Combinatorial Proof of the All Minors Matrix Tree Theorem"
3. **Hodges Determinant:** Hodges - "Eliminating spurious poles from gauge-theoretic amplitudes"
4. **KLT Relations:** Kawai, Lewellen, Tye - "A relation between tree amplitudes..."
5. **Twisted Cohomology:** Mizera - "Scattering Amplitudes from Intersection Theory"

---

*Derivation completed January 2026.*

