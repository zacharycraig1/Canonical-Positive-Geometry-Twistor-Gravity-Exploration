# HODGES FIX APPLIED

## Problem Identified ✅

The issue was that `det(Phi) = 0` for all moment-curve points because:
- **Phi has corank 3** (rank = n - 3 = 3 for n=6)
- Taking a 4×4 determinant of a corank-3 matrix always gives 0
- This is **expected behavior**, not a bug in sampling

## Fix Applied ✅

### 1. Build Full n×n Phi Matrix
- Changed from 4×4 submatrix to full 6×6 matrix
- Off-diagonal: `Phi_{ij} = [i j] / <i j>`
- Diagonal: `Phi_{ii} = - sum_{j != i} Phi_{ij} * (<j x><j y>) / (<i x><i y>)`
  - Uses reference legs x=0, y=5

### 2. Compute Reduced Determinant det'(Phi)
- Remove rows/cols (0,1,2) → keep (3,4,5)
- `det'(Phi) = det(Phi[3,4,5 ; 3,4,5]) / (<01><12><20>)^2`
- This is independent of the choice of removed rows/cols

### 3. Final Amplitude
- `M_6^MHV = det'(Phi) / (∏<i,i+1>)^2`
- Denominator is squared (not just product)

## Removed ❌

- All attempts to compute 4×4 determinants
- Ad-hoc fallback logic to 3×3 minors
- Incorrect assumption that det(Phi) should be nonzero

## Expected Results

1. ✅ `rank(Phi) == 3` (for n=6)
2. ✅ Any 4×4 minor of Phi == 0 (confirms corank-3)
3. ✅ `det'(Phi) != 0` (generic, nonzero)
4. ✅ `det'(Phi)` invariant under choice of removed rows/cols

## Testing

Running `correct_klt_proof.sage` to verify:
- Hodges now returns nonzero values
- KLT matches Hodges (up to normalization)
- Ratio is constant across all points

---

*Status: Fix applied, testing in progress*






