# Pole Accounting for N=6 MHV Gravity

## 1. Physical Singularities
The gravity amplitude $M_6$ should have poles when:
- $s_{ij} \to 0$ (Soft/Collinear)
- $s_{ijk} \to 0$ (Multi-particle factorization)

### 1.1 Soft Limits ($p_s \to 0$)
- For $h=-2$ (negative helicity soft graviton): $M_n \sim \epsilon^{-3}$.
- For $h=+2$ (positive helicity soft graviton): $M_n \sim \epsilon^{+1}$ (vanishes).

### 1.2 Factorization
- On $s_{ijk} \to 0$, $M_6 \to M_L \frac{1}{s_{ijk}} M_R$.
- Residue is product of sub-amplitudes.
- If helicity configurations of $L$ and $R$ are forbidden (e.g. 1 minus in Gravity), residue is zero.

## 2. Origin of Poles in Hodges Formula
The Hodges formula is:
$$ M_n = \text{det}(\Phi_{\text{red}}) \times \langle a b \rangle^8 \times \text{Norm}^{-1} $$
where $\Phi_{ij} = [ij] / \langle ij \rangle$.

### 2.1 Collinear Poles ($<ij> \to 0$)
- If $\langle ij \rangle \to 0$, then $\Phi_{ij}$ and $\Phi_{ji}$ diverge.
- This creates poles in the determinant.
- Since $\det \Phi \sim \Phi_{ij}\Phi_{ji} \sim 1/\langle ij \rangle^2$, this generates the correct gravity scaling $1/p^2$.

### 2.2 Soft Poles
- If $\lambda_s \to \epsilon \lambda_s$:
  - All $\langle i s \rangle \to \epsilon$.
  - All $\Phi_{is} \sim 1/\epsilon$.
  - Diagonal entries $\Phi_{ii}$ also pick up $1/\epsilon$ terms.
- **Issue**: Naive power counting of the determinant suggests $M_6 \sim \epsilon^{-2}$ or worse, even for $h=+2$ where it should be $\epsilon^1$.
- This implies significant cancellations must occur in the determinant for the correct soft behavior to emerge.
- Numerical checks using `hodges_npt_mhv_spinor` show huge values ($1/\epsilon^2$ scaling) for $h=+2$ soft limit, suggesting these cancellations might be numerically unstable or require analytic handling.

## 3. Map to Positive Geometry
- The positive geometry is defined in terms of $z_{ij}$.
- The map is $z_{ij} = \frac{[ij] C_i C_j}{\langle ij \rangle}$.
- Poles in $z_{ij}$ occur when $\langle ij \rangle \to 0$ (Collinear).
- Zeros in $z_{ij}$ occur when $[ij] \to 0$.
- The "Forest Polynomial" $F(z)$ is a polynomial in $z$.
- Singularities of the amplitude arise from:
  1. The map $z \to \text{kinematics}$ (which is singular at collinear limits).
  2. The Jacobian of the pushforward?
  3. The "Canonical Form" structure $\Omega = d \log F$? No, $1/F$.
  - If $F(z)$ has roots, we get poles.
  - But $F(z)$ is a sum of monomials (positive for positive $z$). Roots are at boundaries.

### 3.1 Missing Poles?
- If $F(z)$ is polynomial, where do $1/s_{ijk}$ poles come from?
- They must come from the prefactors or the "integration" measure in the pushforward.
- Current hypothesis: The simple map $z \to p$ might not capture all multi-particle poles manifestly in the geometry of $F$, or they appear as "spurious" singularities of the map that resolve into physical poles.

## 4. Discrepancies Observed
- **Factorization**: Numeric checks confirm $s_{02}$ is a pole, but residue matches only up to a factor (possibly normalization or helicity sum issue).
- **Soft Limit**: Numeric checks fail for $h=+2$ soft limit using Hodges, likely due to cancellation issues described above.




