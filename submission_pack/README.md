# Signed Geometry of MHV Gravity: Forest Expansion and Split Signature

**Main Result:** We identify the sign structure of the forest expansion for $n$-point MHV gravity amplitudes.

The amplitude is computed exactly by the principal minor of a **Weighted Laplacian Matrix** (corank 3), which expands combinatorially into a sum over **3-rooted spanning forests**. We show that this expansion exhibits **signed geometry**—approximately half the forest terms are positive and half are negative, correlating with the $(3,3)$ split signature of the KLT kernel.

**Note:** This work establishes that the natural forest triangulation yields signed (not positive) geometry. Whether an alternative triangulation could yield a positive geometry remains an open question.

---

## 🚀 Reproduction (Verification Suite)

A complete verification suite is available in the **`repro/`** directory.

To run all checks:

```bash
# Linux / Mac
./repro/run_all.sh

# Windows (PowerShell)
.\repro\run_all.ps1
```

### Verification Manifest

**Core Identity & Structure:**
1. **Reference Independence:** `src/scripts/phaseF1_reference_independence.py`
2. **Deletion Independence:** `src/scripts/phaseF2_deletion_set_independence.py`
3. **Pole Order Audit:** `src/scripts/phaseF3_exact_pole_orders.py`
4. **Forest Combinatorics:** `src/scripts/phaseF4_all_minors_forest_expansion.py`
5. **Polytope Geometry:** `src/scripts/phaseF5_newton_polytopes.py`

**Sign Rule & Generalization:**
6. **n=6 Sign Rule:** `src/signed_geometry/verify_chy_sign_derivation.sage`
7. **n=7 Sign Rule:** `src/signed_geometry/generalize_n7.sage`
8. **MTT Consistency:** `tests/test_oracle_match.sage`

---

## 📖 Key Formula

$$ \mathcal{M}_n = (-1)^{n-1} \langle ab \rangle^8 \cdot \frac{\det(\tilde{L}^{(R)})}{\prod_{k \notin R} C_k^2 |\mathcal{N}(R)|^2} $$

Where $\tilde{L}$ is the weighted Laplacian with weights $w_{ij} = [ij]/\langle ij \rangle$ and vertex factors $C_i = \langle i x \rangle \langle i y \rangle$.

**Sign Rule:** Each forest term has sign $\varepsilon(F) = \text{sign}(\prod_e w_e) \times \text{sign}(\prod_v C_v^{\deg(v)})$.

---

## 📂 Structure

- `src/`: Core Python/Sage library (`chy_oracle`, `signed_geometry`).
- `results/`: JSON logs of verification runs.
- `paper/`: Complete LaTeX paper.
- `repro/`: One-click reproduction scripts.
- `verification_packet/`: Full test results with console output.

---

## Key Findings

| n | Forests | Modal Split | Sign Rule Verified |
|---|---------|-------------|--------------------|
| 6 | 108 | (54, 54) | 20/20 samples |
| 7 | 1029 | (515, 514) | 20/20 samples |

The $(54,54)$ and $(515,514)$ splits represent balanced signed geometry—as close to 50/50 as combinatorially possible.
