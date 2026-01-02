# MHV Single Solution Analysis - Execution Checklist

## Pre-Execution Checklist

### Environment Setup
- [ ] SageMath installed (check: `sage --version`)
- [ ] Working directory: `c:\Users\zacha\physics`
- [ ] All required files present:
  - [ ] `src/gravity_proof/analyze_mhv_single_solution.sage` (19 KB)
  - [ ] `src/gravity_proof/amplitude_comparison.sage`
  - [ ] `src/gravity_proof/psi_matrix.sage`
  - [ ] `src/gravity_proof/scattering_solver.sage`
  - [ ] `src/kinematics/spinors.py`
  - [ ] `run_mhv_analysis.ps1`

### Documentation Review
- [ ] Read `RUN_MHV_ANALYSIS_README.md` (quick start)
- [ ] Review `MHV_SINGLE_SOLUTION_INVESTIGATION.md` (technical details)
- [ ] Understand expected outcomes

## Execution

### Run Analysis
```powershell
# Option 1: PowerShell script
.\run_mhv_analysis.ps1

# Option 2: Direct Sage
sage src\gravity_proof\analyze_mhv_single_solution.sage
```

### Monitor Output
Watch for these key sections:
1. ✓ Setup (kinematics, solutions found)
2. ✓ Task 1: Per-solution contributions
3. ✓ Task 2: MHV theory verification
4. ✓ Task 3: Single-solution CHY
5. ✓ Task 4: Pfaffian prefactor check
6. ✓ Task 5: det'(Φ) convention check
7. ✓ Summary and conclusions

## Results Interpretation

### Success Criteria

#### ✅ Primary Success: MHV Theory Confirmed + Exact Match
```
✅ MHV SINGLE SOLUTION THEORY CONFIRMED
   Only 1 of 6 solutions has non-zero Pf'(Ψ)
   ✅ MHV solution matches Hodges amplitude!
   
   CONCLUSION: CHY normalization issue is RESOLVED
```

**Action**: 
- [ ] Document success in proof
- [ ] Update CHY implementation to use only MHV solution
- [ ] Re-run full proof verification
- [ ] Mark issue as RESOLVED

#### ✅ Secondary Success: MHV Theory Confirmed + Normalization Factor
```
✅ MHV SINGLE SOLUTION THEORY CONFIRMED
   Only 1 of 6 solutions has non-zero Pf'(Ψ)
   ⚠️  MHV solution does NOT match Hodges amplitude
   Ratio (single/Hodges) = 1.667
   
✅ Match found with factor 0.6000 (≈ 3/5)
```

**Action**:
- [ ] Document normalization factor: ______
- [ ] Investigate origin of factor (CHY measure? Gauge conventions?)
- [ ] Update `compute_chy_amplitude()` with factor
- [ ] Re-run analysis to verify
- [ ] Mark issue as RESOLVED WITH NORMALIZATION

#### ⚠️ Partial Success: MHV Theory Confirmed, No Factor Match
```
✅ MHV SINGLE SOLUTION THEORY CONFIRMED
   Only 1 of 6 solutions has non-zero Pf'(Ψ)
   ⚠️  MHV solution does NOT match Hodges amplitude
   Ratio (single/Hodges) = 2347.89
   
⚠️  No simple normalization factor found
```

**Action**:
- [ ] Review Pfaffian sign conventions
- [ ] Check det'(Φ) computation
- [ ] Verify Hodges amplitude calculation
- [ ] Consult CHY original papers for conventions
- [ ] Consider alternative gauge fixings

#### ❌ Failure: MHV Theory Not Confirmed
```
⚠️  MHV SINGLE SOLUTION THEORY NOT CONFIRMED
   Found 3 solutions with non-zero Pf'
   (Expected: exactly 1)
```

**Action**:
- [ ] Verify helicity configuration (particles 1,2 neg, 3-6 pos)
- [ ] Check momentum conservation
- [ ] Try different kinematics seed
- [ ] Adjust numerical threshold for zero detection
- [ ] Debug Pfaffian computation

## Post-Execution Analysis

### Data Collection
- [ ] Record dominant solution index: ______
- [ ] Record Pf'(Ψ) at dominant solution: ______
- [ ] Record MHV contribution: ______
- [ ] Record Hodges amplitude: ______
- [ ] Record ratio (MHV/Hodges): ______
- [ ] Record normalization factor (if any): ______

### Verification Steps

If MHV theory confirmed:
- [ ] Verify exactly 5 solutions have |Pf'| < threshold
- [ ] Verify 1 solution has |Pf'| >> threshold
- [ ] Check ratio of 1st to 2nd largest |Pf'| > 10

If normalization factor found:
- [ ] Test factor is rational or simple algebraic number
- [ ] Check if factor relates to (n-3)!, gauge factor, etc.
- [ ] Verify factor is consistent across multiple kinematics

### Multi-Kinematics Test (Robustness)
```python
load("src/gravity_proof/analyze_mhv_single_solution.sage")
all_results = run_multiple_kinematics(num_tests=5, seed_start=0, verbose=False)
```

- [ ] Run with 5 different kinematics
- [ ] Verify MHV theory holds for all
- [ ] Check normalization factor is consistent
- [ ] Record any outliers or failures

## Troubleshooting

### Common Issues

#### Issue: No solutions found
**Symptoms**: "No solutions found!" error
**Fixes**:
- [ ] Try different seed: `main(seed=43)`
- [ ] Check kinematics generation
- [ ] Verify scattering equation setup

#### Issue: Multiple non-zero Pfaffians
**Symptoms**: 2+ solutions with large |Pf'|
**Fixes**:
- [ ] Check if truly MHV (helicities correct?)
- [ ] Adjust threshold: `ZERO_THRESHOLD = 1e-8`
- [ ] Try different kinematics
- [ ] Verify spinor bracket computations

#### Issue: Normalization factor not simple
**Symptoms**: Ratio = 2347.8912... (no simple fraction)
**Fixes**:
- [ ] Check for sign errors (-1 factor?)
- [ ] Verify Pfaffian formula matches CHY convention
- [ ] Check det'(Φ) includes all factors
- [ ] Review gauge-fixing measure

#### Issue: SageMath errors
**Symptoms**: Import errors, syntax errors
**Fixes**:
- [ ] Verify Sage version: `sage --version`
- [ ] Check file paths are correct
- [ ] Ensure all dependencies installed
- [ ] Try loading modules manually in Sage shell

## Documentation

### Update These Files After Success:
- [ ] `PHASE_ATLAS_FIBER_GAUGE_HANDOFF_RESULTS.txt` - Add MHV analysis results
- [ ] Main proof documentation - Note MHV single-solution property
- [ ] `amplitude_comparison.sage` - Add normalization factor if needed
- [ ] Create result summary file

### Record Results:
```markdown
## MHV Single Solution Analysis Results

**Date**: [date]
**Seed**: 42
**MHV Theory**: [Confirmed / Not Confirmed]
**Hodges Match**: [Yes / With Factor / No]
**Normalization Factor**: [value or "None needed"]

### Solution Breakdown:
- Solution 1: Pf' = [value], Contribution = [value]
- Solution 2-6: Pf' ≈ 0

### Conclusion:
[Brief summary of findings and resolution status]
```

## Next Steps Decision Tree

```
Did MHV theory confirm? (1 non-zero solution)
├─ YES
│  ├─ Does MHV solution = Hodges?
│  │  ├─ YES → ✅ RESOLVED
│  │  │  └─ Update CHY implementation
│  │  └─ NO
│  │     └─ Simple normalization factor found?
│  │        ├─ YES → ✅ RESOLVED WITH FACTOR
│  │        │  └─ Apply factor to CHY
│  │        └─ NO → ⚠️ INVESTIGATE CONVENTIONS
│  │           └─ Review CHY papers, check signs
└─ NO → ❌ DEBUG REQUIRED
   └─ Check kinematics, helicities, implementation
```

## Final Checklist

- [ ] Analysis executed successfully
- [ ] Results recorded
- [ ] MHV theory status determined
- [ ] If applicable: Normalization factor identified
- [ ] Documentation updated
- [ ] Next steps identified
- [ ] Implementation updates planned (if needed)
- [ ] Ready to proceed with proof completion

---

**Completion Date**: _______________  
**Analyst**: _______________  
**Status**: [ ] Pass [ ] Pass with Factor [ ] Fail - Needs Debug  
**Notes**: _________________________________________________

