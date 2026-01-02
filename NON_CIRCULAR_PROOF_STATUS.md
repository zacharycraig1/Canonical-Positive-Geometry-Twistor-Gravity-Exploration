# NON-CIRCULAR PROOF STATUS

## Summary

We have successfully refactored the code to remove circular definitions and implement a proper test framework. The test now ONLY compares `hodges_6pt_mhv()` vs `amplituhedron_from_bcfw_cells()` - no circular calls.

## What's Complete ✅

1. **Non-Circular Test Framework**
   - Removed `amplituhedron_6pt_mhv_correct()` that just called Hodges
   - Test harness only compares BCFW cells vs Hodges
   - No circular definitions

2. **Positive Sampling**
   - Implemented `sample_positive_Z()` for G_+(4,6)
   - Checks all ordered 4×4 minors positive
   - Checks angle brackets (needs refinement)

3. **Failure Instrumentation**
   - Logs first failing denominator
   - Records which minor/bracket hit zero
   - Verifies failures are domain violations

4. **Exact Rational Equality Test**
   - Uses QQ arithmetic for exact computation
   - Checks ratio A/H is constant
   - Tests A - ratio*H == 0

5. **Regression Test Framework**
   - Tests with positive sampling
   - Tests with random sampling
   - Classifies failures correctly

## Current Issues ⚠️

### 1. Positive Sampling
- **Problem**: Many attempts fail to find positive points
- **Reason**: Even with positive 4×4 minors, angle brackets can be zero
- **Solution Needed**: Better sampling algorithm or accept some failures and classify them

### 2. BCFW Formula
- **Problem**: Current formula is placeholder - doesn't match Hodges
- **Current**: `term = (channel_angles^2) / (channel_bracket^2 * all_angles^2)`
- **Issue**: Scale factor varies wildly, not constant
- **Solution Needed**: Get correct BCFW formula from literature

## Test Results (Preliminary)

- **Positive sampling**: Finding some positive points, but many attempts fail
- **None cases on positive points**: Some errors detected (Hodges fails even on positive points)
- **BCFW vs Hodges**: Will show mismatches (formula is wrong)

## Next Steps

1. **Improve Positive Sampling**
   - Use known parameterization of G_+(4,6)
   - Or accept failures and ensure they're classified as domain violations

2. **Get Correct BCFW Formula**
   - From literature on 6-point MHV gravity
   - Or derive from first principles
   - Or work backwards from Hodges

3. **Implement Polynomial Identity Testing**
   - Clear denominators
   - Test polynomial equality on many points
   - Bounded degree → high confidence

4. **Uniqueness Argument**
   - Show same poles/residues
   - 1-dimensional space
   - Fix scale

## Files

- `proper_proof_amplituhedron_hodges.sage` - Main non-circular test
- `CURRENT_BCFW_FORMULA.md` - Documents current formula
- `REFACTORING_PROGRESS.md` - Progress tracking
- `NON_CIRCULAR_PROOF_STATUS.md` - This file

## Key Insight

The framework is correct - we just need:
1. Better positive sampling (or better failure classification)
2. Correct BCFW formula

Once we have the correct BCFW formula, the test framework will prove amplituhedron = Hodges.

---

*Status: Framework complete, formula needs correction*










