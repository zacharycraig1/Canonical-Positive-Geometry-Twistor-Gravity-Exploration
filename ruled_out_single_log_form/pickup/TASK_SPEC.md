# Task Specification (Original User Request)

This is the original task specification provided by the user.

---

## CURSOR TASK SPEC — Next Steps (Hodges/KLT → Positive Geometry "Theory" Repo)

**Goal**: Turn the verified Hodges/KLT milestone into a reproducible, test-driven pipeline that can plug into the DCP / positive-geometry constraint-search code.

**Priority order**: Phase 0 → Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 5

**Success criteria**:
- Hodges reduced is invariant under reference spinors and deletion sets on physical kinematics.
- KLT equals Hodges reduced up to the chosen little-group normalization (exact over QQ).
- A clean "oracle" exists that the geometry engine can call.
- A clean mapping exists: spinor-helicity → channel coordinates s_S used by DCP charts.
- Factorization tests exist as a harness to validate candidate forms.

**IMPORTANT**: Work over QQ exactly whenever possible.

---

## Phase 0 — Repo hygiene + Golden tests ✅ COMPLETE

1. `tests/test_hodges_invariance.sage` - Verify Hodges reduced is independent of reference spinors and deletion sets
2. `tests/test_klt_equals_hodges.sage` - Verify exact identity M6^KLT == C(λ) * bar_M6^Hodges
3. `tools/certify_sample.sage` - Certificate generator for reproducibility

---

## Phase 1 — Bridge: spinor-helicity → channel coordinates ✅ COMPLETE

`src/kinematics_map.sage`:
- `spinors_to_sij(lambdas, tilde_lambdas)` → dict[(i,j)] = s_ij
- `sij_to_channels(sij)` → dict[S] = s_S for DCP channels
- `spinors_to_channels(lambdas, tilde_lambdas)` → dict[S] = s_S

---

## Phase 2 — Make Hodges/KLT a callable ORACLE ✅ COMPLETE

`src/oracle_gravity_mhv6.sage`:
- `oracle_M6(lambdas, tilde_lambdas, convention="hodges_reduced"|"klt")`

---

## Phase 3 — Factorization / boundary harness ✅ COMPLETE

1. `src/boundary_sampler.sage`:
   - `sample_near_boundary(S, eps)` → kinematics with s_S ≈ eps

2. `tests/test_factorization_oracle.sage`:
   - Check M6 ~ (M_L * M_R) / s_S near boundaries

---

## Phase 4 — Connect to DCP charts 🔧 INFRASTRUCTURE READY

1. `src/dcp_eval.sage`:
   - `eval_form_on_channels(Omega, channels)` → value
   - `eval_form_on_spinors(Omega, lambdas, tilde_lambdas)` → value

2. `tests/test_dcp_matches_oracle.sage`:
   - Evaluate ratio Ω/oracle across samples

---

## Phase 5 — Uniqueness from consistency ⏳ IN PROGRESS

1. `src/ansatz_builder.sage`:
   - Construct vector space of allowed top log-forms
   - Implement constraints: poles only on channels, projectivity, S6 invariance, factorization

2. `src/solve_constraints.sage`:
   - Solve linear systems exactly over QQ
   - Return surviving basis (should be 1-dimensional)

**Success criterion**: dimension collapses to 1 AND matches oracle on random samples.

---

## Recommended Repo Structure

```
src/
  hodges.sage
  klt.sage
  spinor_sampling.sage
  hodges_sigma.sage
  kinematics_map.sage
  oracle_gravity_mhv6.sage
  boundary_sampler.sage
  dcp_eval.sage
  ansatz_builder.sage
  solve_constraints.sage

tests/
  test_hodges_invariance.sage
  test_klt_equals_hodges.sage
  test_channel_identities.sage
  test_factorization_oracle.sage
  test_dcp_matches_oracle.sage

tools/
  certify_sample.sage
  run_all_tests.sage

results/
  certificates/
  logs/
  figures/
```

