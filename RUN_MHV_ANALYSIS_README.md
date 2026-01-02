# Quick Start: MHV Single Solution Analysis

## Objective
Verify that only ONE of 6 solutions contributes for MHV amplitudes (Du-Teng-Wu theory), resolving the ~2220x CHY vs Hodges discrepancy.

## Quick Run (Recommended)

### Windows PowerShell:
```powershell
.\run_mhv_analysis.ps1
```

### Linux/Mac or Direct Sage:
```bash
sage src/gravity_proof/analyze_mhv_single_solution.sage
```

## What It Does

1. **Solves scattering equations** for 6-point MHV kinematics
2. **Computes Pf'(Ψ)** at each of the 6 solutions
3. **Identifies** which solution(s) have non-zero Pfaffian
4. **Compares** single-solution CHY to Hodges amplitude
5. **Checks** normalization conventions

## Expected Output (Success Case)

```
================================================================================
TASK 1: ANALYZE PER-SOLUTION CONTRIBUTIONS
================================================================================
Solution 1: Pf'(Ψ) = 8.67e+05, Contribution = 7.38e+01 ⭐ DOMINANT
Solution 2: Pf'(Ψ) = 2.34e-12, Contribution = 1.45e-15
Solution 3: Pf'(Ψ) = 7.89e-13, Contribution = 4.23e-16
Solution 4: Pf'(Ψ) = 1.12e-11, Contribution = 8.90e-15
Solution 5: Pf'(Ψ) = 5.67e-12, Contribution = 3.21e-15
Solution 6: Pf'(Ψ) = 9.01e-13, Contribution = 7.11e-16

✅ Single solution dominates (ratio > 10)

================================================================================
TASK 2: VERIFY DU-TENG-WU MHV SINGLE SOLUTION THEORY
================================================================================
Solutions with Pf'(Ψ) ≈ 0: 5
Solutions with Pf'(Ψ) ≠ 0: 1

✅ MHV SINGLE SOLUTION THEORY CONFIRMED!
   Only solution 1 has non-zero Pf'(Ψ)

MHV solution contribution: 7.376e+01
Hodges amplitude:          7.376e+01
Relative difference:       3.45e-09

✅ MHV solution = Hodges (exact match)

================================================================================
SUMMARY AND CONCLUSIONS
================================================================================
✅ MHV SINGLE SOLUTION THEORY CONFIRMED
   Only 1 of 6 solutions has non-zero Pf'(Ψ)
   ✅ MHV solution matches Hodges amplitude!

   CONCLUSION: CHY normalization issue is RESOLVED
```

## If Normalization Factor Needed

If the output shows:
```
Ratio (single/Hodges) = 1.667
```

Look for:
```
Testing normalization factors:
  factor= 0.600: scaled=7.376e+01, rel_diff=5.67e-09 ✅
```

This identifies the missing normalization factor (e.g., 3/5).

## Troubleshooting

### "No solutions found"
- Check kinematics generation (seed might produce degenerate configuration)
- Try different seed: `main(seed=43)`

### "Multiple solutions have non-zero Pf'"
- Verify helicity configuration (should be 1,2 negative; 3,4,5,6 positive)
- Check numerical threshold (may need adjustment)
- This could indicate non-MHV kinematics

### "SageMath not found"
- Install SageMath: https://www.sagemath.org/
- Add to PATH, or use full path to sage executable

## Advanced Usage

### Test Multiple Kinematics
```python
load("src/gravity_proof/analyze_mhv_single_solution.sage")
all_results = run_multiple_kinematics(num_tests=5, seed_start=0, verbose=False)
```

This tests 5 different random kinematics to verify results are generic.

### Customize Analysis
```python
# With specific seed
results = main(seed=123, verbose=True)

# Access results
if results['summary']['mhv_theory_confirmed']:
    print("MHV theory holds!")
```

## Output Files

None by default. Results printed to console.

To save results, modify script to write JSON:
```python
import json
with open('mhv_analysis_results.json', 'w') as f:
    json.dump(results, f, indent=2)
```

## Expected Runtime

- Single kinematics: 1-3 minutes
- Multiple kinematics (5 tests): 5-15 minutes

(Depends on system speed and whether Gröbner basis computation is fast)

## Success Criteria

✅ **MHV Theory Confirmed**: Exactly 1 of 6 solutions has Pf'(Ψ) ≠ 0  
✅ **Hodges Match**: Single solution contribution equals Hodges amplitude (within 1e-8 relative error)  
✅ **Conventions Verified**: Pfaffian and determinant conventions are correct

If all three criteria met → **CHY normalization issue RESOLVED**

## Next Actions After Success

1. Update `compute_chy_amplitude()` to use only MHV solution
2. Apply identified normalization factor (if any)
3. Re-run full proof: `sage src/gravity_proof/run_proof.sage`
4. Verify all success criteria pass

## Questions?

See full documentation: `MHV_SINGLE_SOLUTION_INVESTIGATION.md`

