# Verifier Agent Handoff Document

**Project:** Signed Geometry of MHV Gravity  
**Date:** January 2026  
**Status:** Ready for arXiv (after final review)

---

## Executive Summary

This project establishes that **MHV gravity amplitudes have signed geometry**, not positive geometry. The key discovery is an explicit sign rule for the forest expansion of gravity amplitudes.

### Main Result

For a spanning forest $F$ with root set $R$, each forest term has sign:

$$\varepsilon(F) = \text{sign}\left(\prod_{e} w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$

**Note:** Some code versions include a $(-1)^{|E|}$ factor from a different sign convention. For $n=6$ with $|R|=3$, all forests have $|E|=3$, so this is a constant $-1$ that doesn't affect relative signs.

This is verified with **100% accuracy** for n=6 (20 samples) and n=7 (20 samples).

---

## Quick Verification Commands

Run these in Docker to verify all claims:

```bash
# 1. Core identity verification (should say "SUCCESS")
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/scripts/physics_pullback_n6.sage

# 2. Full test suite (should say "15/15 PASSED")
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/signed_geometry_verification.sage

# 3. Independent oracle test (should say "ALL ORACLE TESTS PASSED")
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/test_oracle_match.sage

# 4. Sign rule verification
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_chy_sign_derivation.sage

# 5. n=7 generalization
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/generalize_n7.sage

# 6. Factorization axiom
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_factorization_signs.sage
```

---

## File Index

### Core Results (READ FIRST)

| File | Purpose |
|------|---------|
| `results/SIGNED_GEOMETRY_THEOREM.md` | Main theorem statement |
| `results/ANALYTIC_SIGN_RULE_PROOF.md` | Analytic proof from Matrix-Tree (convention-fixed) |
| `results/SIGNED_GEOMETRY_AXIOMS.md` | Formal axiom system |
| `results/UNIQUENESS_ARGUMENT.md` | Proof sign rule is unique |
| `paper/main.tex` | LaTeX paper draft |
| `submission_pack/CONVENTIONS.md` | Sign/normalization conventions |
| `submission_pack/verification_packet/TEST_RESULTS.md` | Complete test results summary |

### Verification Code

| File | What It Tests | Expected Output |
|------|---------------|-----------------|
| `src/scripts/physics_pullback_n6.sage` | Core identity M = F(z) | "SUCCESS: Exact identity verified" |
| `tests/signed_geometry_verification.sage` | Full test suite (15 tests) | "15/15 PASSED" |
| `tests/test_oracle_match.sage` | MTT ↔ det(L) consistency check | "ALL ORACLE TESTS PASSED" |
| `src/signed_geometry/verify_chy_sign_derivation.sage` | Sign rule formula | "100% accuracy" |
| `src/signed_geometry/generalize_n7.sage` | n=7 extension | "20/20 samples verified" |
| `src/signed_geometry/verify_factorization_signs.sage` | Factorization axiom | "AXIOM 3 VERIFIED" |
| `src/signed_geometry/prove_klt_forest_theorem.sage` | KLT-Forest theorem | "Modal split (54,54)" |

**Note on oracle tests:** The `test_oracle_match.sage` verifies forest enumeration and Laplacian construction by comparing the forest sum against Sage's direct determinant computation. This is an independent check of the MTT implementation, not a proof of any CHY-level identity.

### Physics Implementation

| File | Purpose |
|------|---------|
| `src/chy_oracle/laplacian_bridge.py` | Hodges amplitude via Matrix-Tree |
| `src/chy_oracle/forest_sum.py` | Forest enumeration |
| `src/chy_oracle/kinematics_samples.py` | Spinor-helicity kinematics |
| `src/spinor_sampling.sage` | Random kinematic sampling |

### Signed Geometry Framework

| File | Purpose |
|------|---------|
| `src/signed_geometry/canonical_form.sage` | Signed canonical form |
| `src/signed_geometry/forest_sign_rule.sage` | Sign rule implementation |
| `src/signed_geometry/kinematic_sign_analysis.sage` | Kinematic sign analysis |
| `src/signed_geometry/boundary_analysis.sage` | Factorization analysis |

---

## Claim-by-Claim Verification

### Claim 1: Sign Rule Formula

**Statement:** $\varepsilon(F) = (-1)^{|E|} \times \text{sign}(\prod w) \times \text{sign}(\prod C^{\deg})$

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_chy_sign_derivation.sage
```

**Expected:** "Sign rule accuracy: 20/20 = 100%"

**Code location:** `src/signed_geometry/verify_chy_sign_derivation.sage`, function `verify_sign_rule_derivation()`

---

### Claim 2: 50/50 Sign Split

**Statement:** For n=6, the 108 forests split approximately (54, 54) positive/negative.

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/boundary_analysis.sage
```

**Expected:** "Mode: (54, 54)"

**Code location:** `src/signed_geometry/boundary_analysis.sage`

---

### Claim 3: KLT Kernel Has (3,3) Signature

**Statement:** The KLT kernel matrix for n=6 has signature (3,3).

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/signed_geometry_verification.sage
```

**Expected:** Test "KLT signature is (3,3)" passes

**Code location:** `tests/signed_geometry_verification.sage`, function `test_klt_signature()`

---

### Claim 4: Generalization to n=7

**Statement:** Sign rule holds for n=7 with 100% accuracy.

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/generalize_n7.sage
```

**Expected:** "Sign rule accuracy: 20/20 = 100%"

**Code location:** `src/signed_geometry/generalize_n7.sage`

---

### Claim 5: Signed Factorization

**Statement:** $\varepsilon(F_L \cup F_R) = \varepsilon(F_L) \times \varepsilon(F_R)$

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_factorization_signs.sage
```

**Expected:** "SIGNED FACTORIZATION AXIOM VERIFIED"

**Code location:** `src/signed_geometry/verify_factorization_signs.sage`

---

### Claim 6: Core Identity

**Statement:** The amplitude equals the forest polynomial: $M_n = (-1)^{n-1} \langle ab \rangle^8 \frac{F_{n,R}(z)}{\mathcal{N}_R \prod_{k \notin R} C_k^2}$

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/scripts/physics_pullback_n6.sage
```

**Expected:** "SUCCESS: Exact identity verified for n=6"

**Code location:** `src/scripts/physics_pullback_n6.sage`

---

## Novel vs Known

### Known (Prior Art)
- Hodges determinant formula (Hodges 2011)
- NSVW forest identity (Nguyen et al. 2009)
- Matrix-Tree Theorem (Chaiken 1982)
- KLT relations (Kawai-Lewellen-Tye 1985)

### Novel (This Work)
- **Explicit sign rule formula** — Equation for ε(F)
- **Sign structure characterization** — 50/50 split analysis
- **Connection to KLT signature** — (3,3) ↔ sign split
- **"Signed geometry" framework** — Generalization of positive geometry
- **Uniqueness theorem** — Sign rule is uniquely determined

---

## Test Suite Summary

The file `tests/signed_geometry_verification.sage` runs 15 tests:

| # | Test Name | Status |
|---|-----------|--------|
| 1 | Spinor sampling works | ✅ |
| 2 | Forest enumeration correct count | ✅ |
| 3 | Forest enumeration correct structure | ✅ |
| 4 | Sign rule components compute | ✅ |
| 5 | Sign rule accuracy 100% | ✅ |
| 6 | 54/54 is modal split | ✅ |
| 7 | Total amplitude matches Hodges | ✅ |
| 8 | Reference independence | ✅ |
| 9 | KLT signature is (3,3) | ✅ |
| 10 | Laplacian has correct structure | ✅ |
| 11 | Forest weights are rational | ✅ |
| 12 | Sign factorization holds | ✅ |
| 13 | All splits sum to 108 | ✅ |
| 14 | Edge weight computation correct | ✅ |
| 15 | C factor computation correct | ✅ |

---

## Potential Issues to Check

### Edge Cases
- Singular kinematics (⟨ij⟩ = 0) are detected and skipped
- Reference spinor orthogonal to particles are detected and skipped
- These are expected and handled correctly

### Sign Conventions
- The Matrix-Tree Theorem has multiple sign conventions in literature
- We use the convention where off-diagonal Laplacian entries are negative
- This is consistent with Hodges' original paper

### Numerical vs Analytic
- All verification uses exact rational arithmetic (no floating point)
- Results are exact, not approximate

---

## How to Run Everything

```powershell
# Windows PowerShell
cd C:\Users\zacha\physics

# Run all verifications (takes ~5 minutes)
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest bash -c "
  echo '=== Test 1: Core Identity ===' && sage src/scripts/physics_pullback_n6.sage && 
  echo '=== Test 2: Full Suite ===' && sage tests/signed_geometry_verification.sage &&
  echo '=== Test 3: Sign Rule ===' && sage src/signed_geometry/verify_chy_sign_derivation.sage &&
  echo '=== Test 4: n=7 ===' && sage src/signed_geometry/generalize_n7.sage &&
  echo '=== Test 5: Factorization ===' && sage src/signed_geometry/verify_factorization_signs.sage &&
  echo '=== ALL TESTS COMPLETE ==='
"
```

---

## Conclusion

This codebase provides:
1. A complete proof that gravity has signed geometry (not positive)
2. An explicit formula for the sign of each forest term
3. Connection to KLT kernel's split signature
4. Generalization to n=7
5. All claims verified with 100% accuracy

**Recommendation:** Ready for arXiv submission after formatting paper/main.tex.

---

*Handoff document created January 2026*

