# Verifier Agent Handoff Document

**Project:** Signed Geometry of MHV Gravity  
**Date:** January 2026  
**Status:** Ready for arXiv submission

---

## Executive Summary

This project establishes that **MHV gravity amplitudes exhibit signed geometry** under the natural forest triangulation. The key discovery is an explicit sign rule for the forest expansion of gravity amplitudes.

**Important:** This work does NOT claim to prove positive geometry for gravity. Rather, it shows the opposite—that the forest triangulation is intrinsically signed, with approximately equal positive and negative contributions.

### Main Result

For a spanning forest $F$ with root set $R$, each forest term has sign:

$$\varepsilon(F) = \text{sign}\left(\prod_{e} w_e\right) \times \text{sign}\left(\prod_v C_v^{\deg(v)}\right)$$

This is verified with **100% accuracy** for n=6 (20 samples) and n=7 (20 samples).

---

## Quick Verification Commands

Run these in Docker to verify all claims:

```bash
# 1. Core identity verification (should say "SUCCESS")
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/scripts/physics_pullback_n6.sage

# 2. Full test suite (should say "15/15 PASSED")
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/signed_geometry_verification.sage

# 3. MTT consistency check (should say "ALL MTT CONSISTENCY TESTS PASSED")
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/test_oracle_match.sage

# 4. Sign rule verification
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_chy_sign_derivation.sage

# 5. n=7 sign rule verification
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
| `results/ANALYTIC_SIGN_RULE_PROOF.md` | Analytic proof from Matrix-Tree |
| `paper/main.tex` | Complete LaTeX paper |
| `submission_pack/CONVENTIONS.md` | Sign/normalization conventions |
| `submission_pack/verification_packet/TEST_RESULTS.md` | Complete test results summary |
| `submission_pack/verification_packet/N7_VERIFICATION.md` | Full n=7 verification logs |

### Verification Code

| File | What It Tests | Expected Output |
|------|---------------|-----------------|
| `src/scripts/physics_pullback_n6.sage` | Core identity M = F(z) | "SUCCESS: Exact identity verified" |
| `tests/signed_geometry_verification.sage` | Full test suite (15 tests) | "15/15 PASSED" |
| `tests/test_oracle_match.sage` | MTT ↔ det(L) consistency | "ALL MTT CONSISTENCY TESTS PASSED" |
| `src/signed_geometry/verify_chy_sign_derivation.sage` | Sign rule formula | "100% accuracy" |
| `src/signed_geometry/generalize_n7.sage` | n=7 extension | "20/20 samples verified" |
| `src/signed_geometry/verify_factorization_signs.sage` | Factorization lemma | "AXIOM 3 VERIFIED" |

---

## Claim-by-Claim Verification

### Claim 1: Sign Rule Formula

**Statement:** $\varepsilon(F) = \text{sign}(\prod w) \times \text{sign}(\prod C^{\deg})$

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_chy_sign_derivation.sage
```

**Expected:** "Sign rule accuracy: 20/20 = 100%"

---

### Claim 2: Balanced Sign Split

**Statement:** For n=6, the 108 forests split approximately (54, 54) positive/negative.

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/boundary_analysis.sage
```

**Expected:** "Mode: (54, 54)"

---

### Claim 3: KLT Kernel Has (3,3) Signature

**Statement:** The KLT kernel matrix for n=6 has signature (3,3).

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/signed_geometry_verification.sage
```

**Expected:** Test "KLT signature is (3,3)" passes

---

### Claim 4: Generalization to n=7

**Statement:** Sign rule holds for n=7 with 100% accuracy.

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/generalize_n7.sage
```

**Expected:** "Sign rule accuracy: 20/20 = 100%", modal split (515, 514)

---

### Claim 5: Signed Factorization

**Statement:** $\varepsilon(F_L \cup F_R) = \varepsilon(F_L) \times \varepsilon(F_R)$ for disjoint subforests

**Verification:**
```bash
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/verify_factorization_signs.sage
```

**Expected:** "SIGNED FACTORIZATION AXIOM VERIFIED"

---

## Novel vs Known

### Known (Prior Art)
- Hodges determinant formula (Hodges 2011)
- NSVW forest identity (Nguyen et al. 2009)
- Matrix-Tree Theorem (Chaiken 1982)
- KLT relations (Kawai-Lewellen-Tye 1985)

### Novel (This Work)
- **Explicit sign rule formula** — Equation for ε(F)
- **Sign structure characterization** — Balanced split analysis
- **Connection to KLT signature** — (3,3) ↔ sign split correlation
- **"Signed geometry" framework** — Generalization of positive geometry concept

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

## Conclusion

This codebase provides:
1. Evidence that gravity exhibits signed geometry (not positive) under the forest triangulation
2. An explicit formula for the sign of each forest term
3. Empirical connection to KLT kernel's split signature
4. Generalization to n=7 with 100% sign rule accuracy
5. All claims verified with exact rational arithmetic

**Note:** Whether an alternative triangulation could yield positive geometry remains an open question.

---

*Handoff document updated January 2026*
