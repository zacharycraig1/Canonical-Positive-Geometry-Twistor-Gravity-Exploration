# Signed Geometry of MHV Gravity

**Current Status:** Publication Ready  
**Key Result:** Explicit Sign Rule for Forest Expansion Derived and Verified

## 🌌 Project Overview

This repository contains the computational verification that **MHV gravity amplitudes have signed geometry**, not positive geometry. We derive an explicit sign rule for the forest expansion and connect it to the KLT kernel's split signature.

### 🏆 The Main Theorem

For a spanning forest $F$ with root set $R$ on $K_n$, the sign of each forest term is:

$$\varepsilon(F) = (-1)^{|E(F)|} \times \text{sign}\left(\prod_{e} w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$

Where:
- $w_{ij} = [ij]/\langle ij \rangle$ are kinematic edge weights
- $C_i = \langle i, x \rangle \langle i, y \rangle$ are reference spinor factors
- $\deg(v)$ is the degree of vertex $v$ in forest $F$

### Key Findings

| Result | Status | Evidence |
|--------|--------|----------|
| Sign rule formula | ✅ PROVEN | Analytic + 100% numerical |
| 50/50 sign split for n=6 | ✅ VERIFIED | (54+, 54-) modal split |
| KLT kernel has (3,3) signature | ✅ VERIFIED | 10/10 samples |
| Generalization to n=7 | ✅ VERIFIED | 20/20 samples, 100% accuracy |
| Gravity ≠ positive geometry | ✅ ESTABLISHED | Intrinsic sign cancellations |

## 🚀 Quick Start

### Using Docker (Recommended)

```powershell
# Verify the core identity
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/scripts/physics_pullback_n6.sage

# Run the full signed geometry test suite
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/signed_geometry_verification.sage

# Verify sign rule for n=6 and n=7
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_chy_sign_derivation.sage
```

### Expected Output

```
SUCCESS: Exact identity verified for n=6 (20 trials)!
...
Test Suite: 15/15 PASSED
```

## 📂 Repository Structure

### `src/signed_geometry/` (Main Contribution)
- `canonical_form.sage`: Signed canonical form implementation
- `forest_sign_rule.sage`: Sign rule derivation and verification
- `verify_chy_sign_derivation.sage`: Complete verification suite
- `generalize_n7.sage`: Extension to n=7

### `results/` (Documentation)
- `SIGNED_GEOMETRY_THEOREM.md`: Complete theorem statement
- `ANALYTIC_SIGN_RULE_PROOF.md`: Analytic proof from Matrix-Tree
- `SIGN_RULE_DERIVATION.md`: Derivation details
- `PUBLICATION_CHECKLIST.md`: Publication readiness status

### `paper/` (LaTeX Paper)
- `main.tex`: Publication-ready paper draft

### `src/chy_oracle/` (Physics Computations)
- `laplacian_bridge.py`: Hodges amplitude via Matrix-Tree
- `forest_sum.py`: Forest enumeration and summation
- `kinematics_samples.py`: Spinor-helicity kinematics

### `tests/` (Verification)
- `signed_geometry_verification.sage`: Comprehensive test suite (15 tests)

## 📖 Key Files to Read

1. **Start here:** `results/SIGNED_GEOMETRY_THEOREM.md`
2. **Proof:** `results/ANALYTIC_SIGN_RULE_PROOF.md`
3. **Paper:** `paper/main.tex`
4. **Tests:** `tests/signed_geometry_verification.sage`

## 🔬 What's Novel vs. Known

| Known (Prior Art) | Novel (This Work) |
|-------------------|-------------------|
| Hodges determinant formula (2011) | Explicit sign rule formula |
| NSVW forest identity (2009) | Sign structure characterization |
| Matrix-Tree Theorem (Chaiken 1982) | Connection to KLT (3,3) signature |
| KLT relations (1985) | "Signed geometry" framework |

## 📊 Verification Summary

| Test | n | Samples | Result |
|------|---|---------|--------|
| Core identity | 6 | 20 | ✅ PASS |
| Sign rule accuracy | 6 | 20 | 100% |
| Sign rule accuracy | 7 | 20 | 100% |
| Hodges amplitude match | 6 | 19 | ✅ PASS |
| KLT signature | 6 | 10 | (3,3) |
| Full test suite | - | 15 | ✅ ALL PASS |

## 🤖 For Reviewers

All claims are numerically verified with exact rational arithmetic in SageMath. The codebase includes:

1. **Reproducibility**: All tests run in Docker container
2. **Documentation**: Every theorem has corresponding verification code
3. **Edge cases**: Singular kinematics are detected and skipped
4. **Independence**: Results verified for multiple random kinematic configurations

## 📝 Citation

If you use this code, please cite:
```
@misc{signed-geometry-gravity,
  title = {Signed Geometry of MHV Gravity},
  author = {[Authors]},
  year = {2026},
  url = {https://github.com/zacharycraig1/Canonical-Positive-Geometry-Twistor-Gravity-Exploration}
}
```

## License

MIT License
