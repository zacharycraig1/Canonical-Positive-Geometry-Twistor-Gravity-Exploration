# Progress: Dimension-1 Candidate Search

## Goal
Find a positive geometry candidate for 6-point MHV gravity amplitudes that:
1. **Factorizes to dimension 1** (unique up to scale) ✓ CURRENT FOCUS
2. **Predicts the correct Hodge structure** (next step)

## Current Configuration

### Script: `54.sage` (dcp_gravity_v15_OPTIMIZED.sage)

**Key Settings:**
- `INTERSECTION_MODE = True` - Using intersection across boundaries
- `INTERSECT_TARGET_DIM = 1` - Target dimension is 1
- `INTERSECT_BOUNDARY_MODE = "ALL_3x3"` - Using all 10 distinct 3|3 boundaries
- `VALIDATION_REQUIRE_DIM1 = True` - Only accept candidates when dim==1
- `INVARIANT_MODE = 'S3xS3'` - Using S3×S3 boundary stabilizer invariants
- `BASE_SEED = 42` - Random seed for reproducibility
- `TOTAL_TRIALS = 100000` - Chart search trials per seed

**Optimizations Applied:**
- PARI kernel (no LLL) - 10-100x speedup
- Parallel workers scaled to all cores
- Modular prefilter with 4 primes
- Sign table precomputation
- Batch matrix operations

## Strategy

The intersection algorithm:
1. Starts with S3×S3 invariant space (typically dim=2-4)
2. Iteratively intersects across 3|3 boundaries
3. Each boundary constraint reduces dimension
4. Stops when dim ≤ 1 or intersection becomes empty
5. Validates robustness with shuffled boundaries and different seeds

## Next Steps

### Phase 1: Find dim=1 candidate (CURRENT)
- Run intersection across all 10 boundaries
- Check if final dimension is exactly 1
- Validate with boundary shuffles and seed offsets
- Save candidate vector `alpha0` when found

### Phase 2: Hodge structure verification (NEXT)
- Extract candidate from dim=1 space
- Compute Hodge numbers/dimensions
- Verify against expected values for 6pt MHV gravity
- Check factorization properties

## Artifacts Saved

When `SAVE_ARTIFACTS = True`, results are saved to:
- `runs/YYYY-MM-DD_HHMMSS_seed42_modeS3xS3/`
  - `config.json` - Configuration used
  - `hit_report_intersection.json` - Full results
  - `candidate_space_basis_in_invariant_coords.sobj` - Transformation matrix T
  - `candidate_space_basis_in_OS_coords.sobj` - Final basis Vinv_final
  - `alpha0_first_basis_vector.sobj` - The dim=1 candidate vector

## Status

- [x] Script configured for dim=1 search
- [x] Intersection mode enabled
- [x] Validation framework in place
- [ ] Dim=1 candidate found
- [ ] Hodge structure verified

## Notes

- The script uses `INTERSECT_STOP_ON_EMPTY = True` - stops if intersection becomes empty
- Can try different invariant modes: 'S6', 'S3xS3', 'S3xS3Z2'
- Boundary order can be randomized for robustness testing
- Multiple seeds can be used for cross-validation












