# Phase D Conventions Ledger

## 1. Integrand Definitions ($M$)

### 1.1 CHY Integrand ($M_{CHY}$)
The standard CHY formula for gravity is:
$$ M_{n} = \int \frac{d\sigma_1 \cdots d\sigma_n}{\text{vol}(SL(2,\mathbb{C}))} {\prod}' \delta(E_a) \det'(\Phi) \det'(\Psi) $$
For MHV, $\det'(\Psi)$ reduces to Parke-Taylor-like factors.

### 1.2 Hodges Determinant Formula ($M_{Hodges}$)
Used in `src/chy_oracle/amplitude_spinor.py`.
$$ M_{Hodges} = \langle 0 1 \rangle^8 \frac{\det(\Phi_{red})}{(\langle r_1 r_2 \rangle \langle r_2 r_3 \rangle \langle r_3 r_1 \rangle)^2} $$
where $\Phi_{red}$ is the reduced Hodges matrix obtained by deleting rows/cols $r_1, r_2, r_3$.
Note: The factor $\langle 0 1 \rangle^8$ accounts for the specific helicity configuration $h_1=h_2=-, h_{3..6}=+$.

## 2. Denominator Spaces & Definitions

### 2.1 Twistor Space ($Z$) - Phase B
- **Space:** Momentum Twistor space $Gr(4,6)$ / $\mathbb{P}^3$.
- **Denominator:**
  $$ D_{\text{found}}(Z) = \frac{(\prod_{i<j}\langle ij\rangle)(\prod_{\text{cyc}}\langle i,i{+}1\rangle)}{\langle 01\rangle^2} $$
  *(Note: The exact form of $D_{found}$ might vary slightly in previous reports, but this captures the structure found).*
- **Numerator:** $N_{\text{found}} = M_{CHY} \cdot D_{\text{found}}$.
- **Status:** Used during the rational function search in Phase B.

### 2.2 Spinor Space ($(\lambda, \tilde\lambda)$) - Phase C3 (Current)
- **Space:** Spinor-Helicity variables (modulo Little Group).
- **Denominator:**
  $$ D_{\text{cyc}}(\lambda) = \left(\prod_{i=1}^n \langle i, i+1 \rangle\right)^2 $$
- **Numerator:** $N_{\text{cyc}} = M_{Hodges} \cdot D_{\text{cyc}}$.
- **Hat Numerator:**
  $$ \hat{N} = \frac{N_{\text{cyc}}}{\langle 0 1 \rangle^8} $$
  This is the object we found to be strictly negative on the positive moment curve.

## 3. Mapping Relationships

The primary relationship to verify is:
$$ M \cdot D_{\text{found}}(Z) \propto (M \cdot D_{\text{cyc}}(\lambda)) \cdot J(\lambda, \tilde\lambda) $$
where $J$ is a Jacobian factor arising from the change of variables or the specific chart choice.

Since $M$ is the physical amplitude (invariant), the relationship is purely between the denominators and the Jacobian.

## 4. Conventions for Phase D
- **Primary Object:** $\hat{N}$ (The "Hat Numerator").
- **Primary Space:** Spinor variables (lifted to Momentum Twistors when needed for geometry).
- **Clearing Factor:** $D_{\text{cyc}}(\lambda)$.
- **Positive Region:** Moment Curve kinematics defined by $Z(t) = (1, t, t^2, t^3)$.

This ledger serves as the single source of truth for Phase D.



