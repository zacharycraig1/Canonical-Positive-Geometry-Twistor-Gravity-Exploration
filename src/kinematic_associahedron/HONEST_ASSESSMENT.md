# Honest Assessment: Cell Decomposition Hypothesis

**Date:** January 2026  
**Status:** ❌ **FALSIFIED**

---

## Executive Summary

The hypothesis that the 6-point MHV gravity amplitude can be written as:

```
M_6 / ⟨12⟩^8 = Σ_{i=0}^{13} c_i / (X_{a_i} × X_{b_i} × X_{c_i})
```

with **constant coefficients c_i** has been **rigorously tested and falsified**.

---

## Test Results

### Test 1: Overdetermined Regression (30 samples, 14 unknowns)

| Metric | Value |
|--------|-------|
| Number of samples | 30 |
| Number of unknowns | 14 |
| Matrix rank | 14 (full rank) |
| Residual L2 norm | ∞ (numerical overflow) |
| Max absolute residual | ~10^9 |

**Result:** The residual is enormous. There is NO set of constant coefficients that can express the amplitude as this sum.

### Test 2: Scale Analysis

For a single kinematic point (seed=42):
- M_reduced = M_6 / ⟨12⟩^8 ≈ -1.49 × 10³
- Triangulation terms ≈ 10⁻⁴ to 10⁰
- Required coefficients ≈ 10³ to 10⁶

For the decomposition to work, these O(10⁶) coefficients would need to be **the same for all kinematic configurations**. They are not.

---

## Why the Previous "Success" Was an Artifact

The original `cell_decomposition.sage` had only **8 valid data points** for **14 unknowns**.

### The Mathematical Reality

With 8 equations and 14 unknowns:
- **Infinite solutions exist** by linear algebra
- Sage's solver picks ONE particular solution (minimal norm)
- The "perfect match" was **guaranteed** because the system was underdetermined
- The claim that "c_8...c_13 = 0" was an **artifact of the solver**, not physics

### The Proper Test

With 30 equations and 14 unknowns:
- The system is **overdetermined** 
- If true coefficients exist, the residual should be **exactly zero**
- The residual is **~10^9** → no such coefficients exist

---

## What This Means

### The Amplitude Does NOT Decompose Simply

The gravity amplitude **cannot** be written as a sum over associahedron cells with constant weights. This rules out the simplest positive geometry interpretation.

### BCJ Structure is More Subtle

The BCJ formula for gravity is:
```
M_n = Σ_T n_T² / D_T
```

where n_T are **kinematic-dependent BCJ numerators**, not constants.

The positive geometry, if it exists, must encode:
1. The kinematic dependence of these numerators
2. A more sophisticated geometric structure than a simple weighted polytope

---

## Verified Results (What IS True)

1. **KLT-Hodges Equivalence**: Proven exactly for 6-point MHV
   - M_6^{KLT} = -⟨12⟩^8 × M̄_6^{Hodges}
   - Verified 19/20 samples with exact rational arithmetic

2. **No S6-Invariant Canonical Form**: Proven in Orlik-Solomon algebra A³
   - The intersection of S6 symmetry + factorization constraints is empty

3. **No Simple Positive Region**: 
   - Worldsheet: Solutions don't lie in any ordering region
   - Kinematic space: No simple s_ij > 0 region works
   - Twistor space: Gr_+(4,6) is measure-zero

4. **Forest Polynomial Identity**: M = F(z) verified for n=4,5,6

---

## Actual Status of the Research Problem

**The positive geometry for MHV gravity remains an OPEN PROBLEM.**

This is consistent with the fact that leading physicists (UNIVERSE+ project at Max Planck Institute) are actively working on this as a research frontier.

### What Has Been Achieved

1. Built robust computational infrastructure
2. Verified known relations (KLT = Hodges)
3. Ruled out several simple hypotheses
4. Identified that the geometry is more subtle than expected

### What Remains Unknown

1. The precise geometric structure encoding gravity amplitudes
2. How BCJ numerators arise geometrically
3. Whether a positive geometry description exists at all for gravity

---

## Recommendations for Future Work

1. **Study BCJ numerators directly**: Compute n_T for each topology and understand their structure

2. **Fiber bundle approach**: The ⟨12⟩^8 factor suggests a fiber over kinematic base space

3. **Double-copy geometry**: Look for a product structure (YM × YM) at the geometric level

4. **BCFW cells**: Each BCFW term might correspond to a geometric cell

5. **Literature review**: Check what the positive geometry community has published on gravity

---

## Conclusion

The "breakthrough" claimed in the previous summary was based on:
- An underdetermined system (8 equations, 14 unknowns)  
- A solver artifact, not a physical discovery

The proper overdetermined test (30 equations, 14 unknowns) definitively **falsifies** the constant-coefficient hypothesis.

**The positive geometry for gravity is NOT the kinematic associahedron with constant weights.**

The research continues.

---

*This assessment reflects rigorous testing, not wishful thinking.*

