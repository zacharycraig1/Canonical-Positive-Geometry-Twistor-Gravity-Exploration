# STATUS: Correct KLT Implementation

## Current Issue

**Problem**: Hodges formula returns 0 for all moment-curve points
- `det(Phi) = 0` for moment-curve configurations
- This makes the entire amplitude 0

**Possible Causes**:
1. Moment-curve points are in a special configuration that makes Phi singular
2. The Hodges formula needs a different normalization
3. The Phi matrix construction has a bug

## Fixes Applied

1. ✅ Fixed `simplify_rational()` bug (Rational objects don't have this method)
2. ✅ Added fallback to 3×3 minors if det(Phi) = 0
3. ⏳ Testing if this resolves the issue

## Next Steps

1. **If det(Phi) still = 0**: Investigate why moment-curve points make Phi singular
   - Check if this is expected behavior
   - Try different particle choices for Phi matrix
   - Check normalization factors

2. **If det(Phi) ≠ 0**: Verify KLT matches Hodges
   - Check ratio is constant
   - Verify exact equality

3. **Alternative**: Use different positive sampling method
   - Not moment-curve, but still guaranteed positive
   - See if Phi is non-singular for those points

## Implementation Status

- ✅ Moment-curve sampler
- ✅ Tight domain checks
- ✅ Correct KLT momentum kernel
- ✅ Correct KLT gravity amplitude
- ✅ Exact equality test
- ⏳ Hodges formula (det(Phi) = 0 issue)

---

*Status: Testing fix for det(Phi) = 0*










