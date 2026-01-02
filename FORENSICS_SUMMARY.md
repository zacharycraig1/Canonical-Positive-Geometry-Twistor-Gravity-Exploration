# FORENSICS SUMMARY: Understanding the "Failures"

## Key Finding

**All failures are due to degenerate kinematics, NOT formula mismatches.**

## What We Discovered

### 1. The "Direct Method" is Circular
In `correct_amplituhedron_hodges.sage`:
```python
def amplituhedron_6pt_mhv_correct(twistor):
    # ...
    return hodges_6pt_mhv(twistor)  # They're literally the same function!
```

This explains why we get 193/200 matches - they're computing the exact same thing. This is not a proof, it's a definition.

### 2. The BCFW Method Needs Work
The BCFW cell sum does NOT equal Hodges with the current formula:
- Only 7 matches (with scaling)
- Large errors (min=0.99, avg=3333)
- Scale factor varies wildly

**The BCFW formula is incorrect and needs to be fixed.**

### 3. All Failures Are Degenerate Cases
Forensics showed:
- 28 out of 200 points returned `None`
- 100% of failures due to:
  - Zero angle brackets: `zero_angle_bracket`
  - Zero square brackets: `zero_square_bracket`
- All failing points have `positivity_ok = False`
- These are not in the amplituhedron domain

## What This Means

### The Good News ✅
- When both A and H are defined (non-degenerate), they match exactly
- No true formula mismatches
- Failures are expected (degenerate kinematics)

### The Bad News ❌
- The "direct method" is not a proof (it's circular)
- The BCFW method doesn't work (formula is wrong)
- We need the correct BCFW formula for gravity MHV

## What We Need to Do

### 1. Fix the BCFW Formula
The current BCFW cell contribution is:
```python
term = (channel_angles^2) / (channel_bracket^2 * all_angles^2)
```

This is wrong. We need the correct BCFW formula for 6-point MHV gravity from literature.

### 2. Sample from Positive Region
Instead of random integers, sample from:
- Positive 4×4 minors: `⟨i j k l⟩ > 0` for `i < j < k < l`
- Positive angle brackets: `⟨i i+1⟩ > 0`

This will eliminate the degenerate cases.

### 3. Prove from First Principles
Once we have the correct BCFW formula:
- Show BCFW sum = Hodges on positive generic points
- Use uniqueness argument (same poles, same residues, 1-dimensional)
- Or polynomial identity testing (clear denominators, check equality)

## The Real Question

**Does the amplituhedron (computed from BCFW cells) equal Hodges?**

Current answer: **We don't know yet** because:
1. The BCFW formula is wrong
2. We need to fix it first
3. Then test on positive points

## Next Steps

1. **Find correct BCFW formula** for 6-point MHV gravity
2. **Implement positive sampling** to avoid degeneracies
3. **Test BCFW sum = Hodges** on positive points
4. **Prove rigorously** using uniqueness or polynomial identity testing

## Conclusion

The forensics revealed:
- ✅ Failures are all degenerate (expected)
- ❌ BCFW formula is wrong (needs fixing)
- ❌ "Direct method" is circular (not a proof)

**We need to fix the BCFW formula and test properly on positive points to get a real proof.**

---

*Status: Forensics complete. BCFW formula needs correction.*










