#!/usr/bin/env sage
# Optimized Boundary Search - Physics Breakthrough Strategy
# Uses 4 well-chosen, independent boundaries as recommended by mathematical insights
# Theory: Start with fewer boundaries, add more if dim > 1
# This avoids over-constraining and finding empty intersections

load('54.sage')

# KEY OPTIMIZATION: Use 4 strategic, independent boundaries
# These cover different aspects of the geometry and avoid redundancy
INTERSECT_BOUNDARY_MODE = "CUSTOM"
INTERSECT_BOUNDARIES = [
    (1,2,3),  # Standard boundary - splits {1,2,3}|{4,5,6}
    (1,2,4),  # Different perspective - covers different channel combinations
    (1,3,5),  # Maximally independent - spreads across points
    (2,4,6),  # Complementary structure - covers remaining aspects
]

# These boundaries are chosen to be:
# - Independent (not redundant)
# - Cover different geometric aspects
# - Balance left/right splits
# - Avoid over-constraining (which could lead to empty intersection)

# Keep multi-strategy search enabled
MULTI_STRATEGY_SEARCH = True
STRATEGY_INVARIANT_MODES = ['S3xS3', 'S3xS3Z2', 'S6']  # Try all modes
STRATEGY_STOP_ON_DIM1 = True  # Stop when dim=1 found

# Optimized search parameters
TOTAL_TRIALS = 50000  # Balanced - fast but thorough
MAX_SOLUTIONS = 100   # Enough to find good charts

# Validation - keep enabled for physics breakthrough verification
RUN_CANDIDATE_VALIDATION = True
VALIDATION_SHUFFLES = 1
VALIDATION_SEED_OFFSETS = [0]
VALIDATION_REQUIRE_DIM1 = True  # Only accept dim=1 candidates

log("\n" + "="*70)
log("OPTIMIZED BOUNDARY SEARCH - Physics Breakthrough Strategy")
log("="*70)
log(f"Using {len(INTERSECT_BOUNDARIES)} strategic boundaries:")
for i, b in enumerate(INTERSECT_BOUNDARIES, 1):
    log(f"  {i}. {b}")
log("="*70)

if __name__ == '__main__':
    main()



