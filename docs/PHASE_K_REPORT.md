# Phase K Report: Physical Factorization & Geometric Pushforward

## 1. Executive Summary
Phase K focused on bridging the gap between the formal geometric object (Forest Polytope) and physical gravity amplitudes. We implemented precise checks for physical factorization, soft/collinear limits, and the 'Stringy Canonical Form' pushforward.

**Key Achievements:**
- **Physical Factorization (K1):** Implemented src/scripts/test_factorization_physical_n6.sage to verify the residue of M6 on the s123 channel matches ML x MR.
- **Field Theory Limits (K2):** Implemented src/scripts/test_soft_limit_n6.sage (Soft Graviton) and src/scripts/test_collinear_limit_n6.sage (Collinear Splitting).
- **Geometric Pushforward (K3):** Implemented Route K3-A via src/posgeom/stringy_form.py, computing the canonical form of the Forest Polytope using the saddle-point method on the 'Stringy' integral.
- **Documentation (K4):** Clarified the novelty of the 'Forest Polytope Pushforward' vs standard Matrix-Tree theorems.

## 2. Detailed Results

### 2.1 Physical Factorization (K1)
*Script:* src/scripts/test_factorization_physical_n6.sage
- **Method:** Deforms momentum twistors to enforce < 5 0 2 3 > -> 0 (implying s012 -> 0).
- **Verification:** Computes Res = lim s012 M6 and compares with M4(0,1,2,-P) x M4(P,3,4,5).
- **Status:** Implemented. Run script to generate RESULTS/factorization_n6.txt.

### 2.2 Soft & Collinear Limits (K2)
*Scripts:* src/scripts/test_soft_limit_n6.sage, src/scripts/test_collinear_limit_n6.sage
- **Soft Limit:** Checks M6 -> S(0) M5 as p5 -> 0. Uses Weinberg soft factor.
- **Collinear Limit:** Checks scaling of M6 as 4 || 5.
- **Status:** Implemented. Run scripts to verify universality.

### 2.3 Forest Polytope Pushforward (K3)
*Code:* src/posgeom/stringy_form.py, src/scripts/test_pushforward_forest.sage
- **Theory:** The Amplitude is identified as the canonical form of the Newton Polytope of the 'Forest Polynomial' F(X) = sum (prod z_e) X^F.
- **Implementation:**
  - StringyCanonicalForm class solves saddle point equations grad log F(X) ~ s/X.
  - Computes contribution 1 / det(H_eff).
  - Matches against canonical_polytope (geometric volume of dual).
- **Conjecture:** The pushforward of the canonical form of the Forest Polytope (under the scattering equations map) is exactly the MHV Gravity Amplitude.

## 3. Main Theorem Candidate

**Theorem (Conjecture):**
The n-point MHV Gravity Amplitude M_n, MHV is the **Pushforward** of the Canonical Form of the **3-Rooted Spanning Forest Polytope** P_n,R to the space of kinematic variables z_ij, via the scattering equations associated with the potential W = sum (prod z_ij) prod x_ij.

**Correspondence:**
| Geometry | Physics |
|---|---|
| Object: Newton Polytope of Forest Polynomial | Scattering Amplitude |
| Facets of Polytope | Factorization Channels |
| Canonical Form Omega(P) | Amplitude Value |
| Variables z_ij (Coefficients) | Kinematics ([ij]/< ij >) Ci Cj |

## 4. Next Steps
- Run the provided scripts to populate RESULTS/ tables.
- Refine the regularization parameter s in the Stringy Form calculation to extract the exact rational number.
- Extend to N=7 to verify scaling of complexity.
