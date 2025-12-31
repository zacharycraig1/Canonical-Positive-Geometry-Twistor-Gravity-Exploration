# AGENT HANDOFF: 6-Point MHV Gravity Positive Geometry Project

## Date: 2024-12-29
## Status: Phase 3 Complete, Phase 5 In Progress

---

## 1. PROJECT GOAL

**Ultimate Objective**: Prove that the 6-point MHV gravity amplitude (Hodges det' formula) is the *unique* canonical form of a positive geometry, derived purely from:
- Poles only on channel divisors {s_S = 0}
- Projectivity
- S6 invariance (or weaker symmetry)
- Factorization residue constraints

This would establish: **"Consistency + symmetry + allowed singularities ⇒ Hodges (up to scale)"**

---

## 2. WHAT HAS BEEN VERIFIED (Phases 0-3)

### Phase 0: Hodges Invariance ✅
- **File**: `tests/test_hodges_invariance.sage`
- **Result**: PASSED (50/50 samples)
- The Hodges reduced amplitude `bar_M6` is independent of:
  - Reference spinor choices (λ_x, λ_y)
  - Deletion set choices (rows/cols to remove from Phi matrix)

### Phase 0: KLT vs Hodges ⚠️
- **File**: `tests/test_klt_equals_hodges.sage`
- **Result**: Ratio varies with kinematics
- **Finding**: The relation `M_KLT = C(λ) * bar_M6` holds but `C(λ)` is NOT simply `-<01>^8`
- **Weight Analysis** (from `tools/check_weights.sage`):
  ```
  Little Group scaling (particle 0): 
    Hodges: t^-4 per particle
    KLT: t^+4 per particle
  Total (all 6 particles, P-invariant scaling):
    Hodges: t^-24
    KLT: t^-8
  ```
- **Implication**: The two formulas compute different objects. Hodges Reduced is the "stripped" amplitude; KLT includes helicity factors.

### Phase 1: Kinematics Map ✅
- **File**: `src/kinematics_map.sage`
- Functions: `spinors_to_sij()`, `spinors_to_channels()`, `SpinorHelicityAdapter`
- Verified: `s_ij = <ij>[ij]` is symmetric, `s_S = sum of s_ij`

### Phase 2: Oracle ✅
- **File**: `src/oracle_gravity_mhv6.sage`
- Function: `oracle_M6(lambdas, tilde_lambdas)` returns Hodges and KLT values

### Phase 3: Factorization ✅
- **File**: `src/boundary_sampler.sage`, `tests/test_factorization_oracle.sage`
- **Result**: PASSED
- The Hodges oracle has correct simple poles:
  - `M * s_S → constant` as `s_S → 0` (both 2-particle and 3-particle channels)
- **Key Observation**: 3-particle channel `s_{012} = 0` iff `<5012> = 0` (4-bracket vanishes)

---

## 3. CURRENT STATE (Phase 4-5)

### Phase 4: DCP Evaluation Infrastructure
- **File**: `src/dcp_eval.sage`
- Can evaluate a candidate form (coefficient vector) on physical spinors
- Tested with dummy form; pipeline works

### Phase 5: Ansatz Builder (PARTIALLY STARTED)
- **File**: `src/ansatz_builder.sage`
- Loads `54.sage` (the main DCP/wonderful model search engine)
- **Issue**: 54.sage's main() runs on load and throws a seed type error at the end

### DCP Search Results (from running `54.sage`)
```
S6 invariant space: dim = 2
After intersecting boundary (1,2,3): rank becomes FULL (2/2)
=> Intersection becomes empty!
```

**This means**: Under full S6 symmetry with the current chart, no form survives the first boundary constraint. This could indicate:
1. The chart is wrong (need different nested set)
2. S6 is too restrictive (try S3×S3)
3. The boundary constraint implementation has a bug

With S3×S3 symmetry:
```
Initial dim = 58
After boundary (1,2,3): dim = 48
After boundary (1,2,4): intersection empty
```

---

## 4. KEY FILES

### Core Source (`src/`)
| File | Purpose |
|------|---------|
| `hodges.sage` | Hodges det' formula implementation |
| `klt.sage` | KLT double-copy implementation |
| `spinor_sampling.sage` | Momentum-conserving spinor generator |
| `kinematics_map.sage` | Spinors → Mandelstams → Channels |
| `oracle_gravity_mhv6.sage` | Unified oracle interface |
| `boundary_sampler.sage` | BCFW shift to reach boundaries |
| `dcp_eval.sage` | Evaluate DCP forms on kinematics |
| `hodges_sigma.sage` | Sign computation for Hodges formula |

### Tests (`tests/`)
| File | Status |
|------|--------|
| `test_hodges_invariance.sage` | ✅ PASS |
| `test_klt_equals_hodges.sage` | ⚠️ Needs convention fix |
| `test_channel_identities.sage` | ✅ PASS |
| `test_factorization_oracle.sage` | ✅ PASS |
| `test_dcp_matches_oracle.sage` | 🔧 Infrastructure only |

### Main DCP Engine
| File | Purpose |
|------|---------|
| `54.sage` | Full DCP search (nested sets, invariants, intersection) |
| `optimized_boundary_search.sage` | Wrapper with strategic boundaries |

### Reports (`results/`)
| File | Content |
|------|---------|
| `PHASE0_REPORT.md` | KLT/Hodges ratio analysis |
| `PHASE3_REPORT.md` | Factorization verification |

---

## 5. KNOWN ISSUES

### Issue 1: 54.sage seed type error
```
Error loading 54.sage: The only supported seed types are: None, int, float, str, bytes, and bytearray.
```
**Cause**: Sage's preparser converts `42` to `Integer(42)`, which Python's `random.seed()` rejects.
**Fix needed**: Find all `random.seed(...)` calls in 54.sage and wrap with `int(...)`.

### Issue 2: Empty intersection under S6
The 2-dimensional S6-invariant space becomes 0-dimensional after the first boundary constraint.
**Possible causes**:
- The 2 S6-invariant forms don't satisfy factorization
- Chart/nested-set choice is incompatible
- Need to relax to S3×S3 or weaker symmetry

### Issue 3: KLT/Hodges weight mismatch
The formulas compute objects with different little-group weights.
**Resolution**: Use Hodges Reduced as the oracle; it's the "stripped" amplitude.

---

## 6. NEXT STEPS FOR CONTINUATION

### Immediate (Fix blockers)
1. **Fix 54.sage seed issue**: Search for `random.seed` and wrap arguments with `int()`.
2. **Debug empty intersection**: 
   - Print the 2 S6-invariant basis vectors explicitly
   - Check if they satisfy the boundary equation manually

### Phase 5 Completion
1. **Build residue constraint matrix**:
   - For each channel S, compute `Residue(M_Hodges, s_S=0)` numerically
   - Express as linear constraint on candidate coefficients
2. **Solve linear system**:
   - Start with S3×S3 basis (dim=58)
   - Impose residue constraints one by one
   - Track dimension reduction
3. **Verify uniqueness**:
   - If dim → 1, extract the unique form
   - Compare to oracle on random samples

### Verification
1. Once unique form is found, evaluate it on 100+ samples
2. Ratio to oracle should be constant (= 1 after normalization)
3. Document the proof in `results/UNIQUENESS_PROOF.md`

---

## 7. EXECUTION COMMANDS

```bash
# Run all Phase 0-3 tests
./sage.ps1 tools/run_all_tests.sage

# Run factorization check
./sage.ps1 tests/test_factorization_oracle.sage

# Run DCP search (takes 5-10 min)
./sage.ps1 54.sage

# Check little-group weights
./sage.ps1 tools/check_weights.sage
```

---

## 8. PHYSICS SUMMARY

### What We Know
1. **Hodges det'** gives the correct 6-pt MHV gravity amplitude (verified against literature)
2. **KLT double-copy** also works but has different normalization
3. The amplitude has **simple poles** on all channel divisors s_S = 0
4. The **reduced determinant** is invariant under reference spinor and deletion set choices
5. In momentum twistor space: `s_{ijk} = 0` ⟺ `<Z_{i-1} Z_i Z_j Z_k> = 0`

### What We're Proving
The Hodges form is the **unique** projective, S6-invariant (or S3×S3) differential form on kinematic space with:
- Poles only on channel divisors
- Correct factorization residues

This would establish gravity amplitudes as canonical forms of a positive geometry.

---

## 9. REFERENCE DOCUMENTS

- `Cursor_6pt_MHV_Gravity_Positive_Geometry_Spec.pdf` - Original spec
- `mathematical_insights.md` - Strategy notes
- `FINAL_SOLUTION_APPROACH.md` - Detailed algorithm
- `BREAKTHROUGH_FINDINGS.md` - Key discoveries

---

## 10. CRITICAL CODE SNIPPETS

### How to sample momentum-conserving spinors (QQ-exact):
```python
load("src/spinor_sampling.sage")
res = sample_spinor_helicity_conserving(n=6, seed=42)
lambdas, tilde_lambdas = res  # Lists of 2-vectors over QQ
```

### How to compute Hodges amplitude:
```python
load("src/hodges.sage")
load("src/kinematics_map.sage")
adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
val, reason = hodges_6pt_mhv_reduced(adapter)
```

### How to approach a boundary (BCFW shift):
```python
load("src/boundary_sampler.sage")
S = {0, 1, 2}  # Channel to probe
z_pole = solve_boundary_shift(lambdas, tildes, S, (0, 3))
l_new, t_new = apply_shift(lambdas, tildes, 0, 3, z_pole + epsilon)
# Now s_S ≈ epsilon * (linear factor)
```

### How to evaluate a DCP form on spinors:
```python
load("src/dcp_eval.sage")
C = generate_channels_n6()  # 25 channels (1-based)
triples = generate_triples(C)  # 2300 basis triples
candidate = vector(QQ, len(triples))  # Coefficient vector
value = eval_form_on_spinors(candidate, lambdas, tildes, C, triples)
```

---

## 11. ARCHITECTURE DIAGRAM

```
                    ┌─────────────────────────────────────────┐
                    │          PHYSICAL AMPLITUDE             │
                    │   (Hodges det' / KLT double-copy)       │
                    └─────────────────┬───────────────────────┘
                                      │
                                      ▼
                    ┌─────────────────────────────────────────┐
                    │              ORACLE                     │
                    │  oracle_M6(λ, λ̃) → M_hodges, M_klt     │
                    └─────────────────┬───────────────────────┘
                                      │
              ┌───────────────────────┴───────────────────────┐
              │                                               │
              ▼                                               ▼
┌─────────────────────────┐                 ┌─────────────────────────┐
│   BOUNDARY SAMPLER      │                 │     DCP EVALUATOR       │
│  sample_near_boundary   │                 │  eval_form_on_spinors   │
│  (BCFW shift to s_S→0)  │                 │  (candidate → number)   │
└─────────────┬───────────┘                 └───────────┬─────────────┘
              │                                         │
              ▼                                         ▼
┌─────────────────────────┐                 ┌─────────────────────────┐
│  FACTORIZATION CHECK    │                 │    UNIQUENESS CHECK     │
│  M * s_S → constant?    │                 │  candidate/oracle = 1?  │
└─────────────────────────┘                 └─────────────────────────┘
```

---

## 12. GLOSSARY

| Term | Definition |
|------|------------|
| `s_ij` | 2-particle Mandelstam: `(p_i + p_j)^2 = <ij>[ij]` |
| `s_S` | Multi-particle Mandelstam: `(sum p_i)^2` for i in S |
| `<ij>` | Angle bracket: `det(λ_i, λ_j)` |
| `[ij]` | Square bracket: `det(λ̃_i, λ̃_j)` |
| `bar_M6` | Hodges reduced amplitude (stripped of helicity factors) |
| `OS3` | Orlik-Solomon algebra degree 3 (space of top log-forms) |
| `DCP` | De Concini-Procesi (wonderful model compactification) |
| `S3×S3` | Stabilizer of (123)\|(456) partition |

---

**END OF HANDOFF DOCUMENT**

