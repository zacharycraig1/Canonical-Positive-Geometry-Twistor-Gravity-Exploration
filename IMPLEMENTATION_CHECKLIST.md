# Implementation Checklist

## Phase 1: Fix Pushforward (Week 1)

### Setup
- [x] Create `src/debug/` directory
- [x] Create `src/debug/__init__.py`
- [x] Create `src/debug/pushforward_n6_analysis.sage`

### Diagnostic Tasks
- [ ] Run `pushforward_n6_analysis.sage` to test ratio constancy
- [ ] Analyze Jacobian eigenvalue structure
- [ ] Compare n=4, n=5, n=6 behavior
- [ ] Document findings in `PHASE1_REPORT.md`

### Decision Point
- [ ] Determine if ratio is constant (→ find normalization)
- [ ] Determine if ratio varies (→ proceed to Phase 2)

---

## Phase 2: BCFW Amplituhedron (Weeks 2-3)

### Setup
- [x] Create `src/amplituhedron/` directory
- [x] Create `src/amplituhedron/__init__.py`
- [x] Create `src/amplituhedron/momentum_twistor.sage`
- [x] Create `src/amplituhedron/bcfw_cells.sage`
- [x] Create `src/amplituhedron/verify_bcfw_hodges.sage`
- [x] Create `src/amplituhedron/run_verification.sage`

### Core Implementation
- [ ] Verify momentum twistor brackets are correct
- [ ] Verify BCFW cell enumeration (Catalan numbers)
- [ ] Implement correct cell canonical form formula
- [ ] Test against Hodges determinant

### Verification
- [ ] Run `sage src/amplituhedron/run_verification.sage`
- [ ] Analyze BCFW vs Hodges discrepancy
- [ ] Identify missing normalization factors
- [ ] Document results in `PHASE2_REPORT.md`

### Success Criteria
- [ ] BCFW sum = Hodges for n=6 (exact or up to constant)
- [ ] Positivity region verified
- [ ] Formula documented mathematically

---

## Phase 3: Intersection Theory (Week 4, if needed)

### Setup
- [x] Create `src/intersection/` directory
- [ ] Create `src/intersection/scattering_eqs.sage`
- [ ] Create `src/intersection/chy_gravity.sage`

### Implementation
- [ ] Implement scattering equation solver
- [ ] Implement CHY gravity integrand
- [ ] Connect to forest polynomial
- [ ] Document in `PHASE3_REPORT.md`

---

## Code Quality Checklist

### Every Session
- [ ] Run smoke tests before starting
- [ ] Use exact arithmetic (QQ) for verification
- [ ] Commit changes to version control
- [ ] Update this checklist

### Every Function
- [ ] Document mathematical formula in comments
- [ ] Add unit test
- [ ] Handle edge cases (singular kinematics)

---

## Key Commands

```bash
# Smoke test
sage tests/test_n4_n5.sage

# Phase 1 diagnostic
sage src/debug/pushforward_n6_analysis.sage

# Phase 2 verification
sage src/amplituhedron/run_verification.sage

# Full test suite
python run_core_tests.py
```

---

## Files Created

### Phase 1
- `src/debug/__init__.py`
- `src/debug/pushforward_n6_analysis.sage`

### Phase 2
- `src/amplituhedron/__init__.py`
- `src/amplituhedron/momentum_twistor.sage`
- `src/amplituhedron/bcfw_cells.sage`
- `src/amplituhedron/verify_bcfw_hodges.sage`
- `src/amplituhedron/run_verification.sage`

### Documentation
- `GRAVITY_AMPLITUHEDRON_PLAN.md`
- `IMPLEMENTATION_CHECKLIST.md` (this file)

---

*Last updated: December 31, 2025*


