#!/usr/bin/env sage
# Targeted search: Start with fewer boundaries, add more if needed
# More efficient approach - finds dim=1 faster

load('54.sage')

# OPTIMIZATION: Start with strategic boundary selection
# Theory: Some boundaries are more constraining than others
# Start with a minimal set that should reduce to dim=1

# Strategy 1: Start with 4 well-chosen boundaries
INTERSECT_BOUNDARY_MODE = "CUSTOM"
INTERSECT_BOUNDARIES = [
    (1,2,3),  # Standard boundary
    (1,2,4),  # Different split
    (1,3,5),  # Another perspective
    (2,4,6),  # Complementary
]

# If this doesn't reach dim=1, we'll add more
TOTAL_TRIALS = 50000
MAX_SOLUTIONS = 100
MULTI_STRATEGY_SEARCH = True
STRATEGY_STOP_ON_DIM1 = True

# Faster validation
VALIDATION_SHUFFLES = 0
VALIDATION_SEED_OFFSETS = []
RUN_CANDIDATE_VALIDATION = False

def adaptive_boundary_selection():
    """Adaptively add boundaries if dim > 1."""
    global INTERSECT_BOUNDARIES
    
    # All possible boundaries
    all_boundaries = all_distinct_3x3_boundaries(6)
    
    # Start with 4
    if len(INTERSECT_BOUNDARIES) < 4:
        INTERSECT_BOUNDARIES = all_boundaries[:4]
    
    log(f"Starting with {len(INTERSECT_BOUNDARIES)} boundaries: {INTERSECT_BOUNDARIES}")

if __name__ == '__main__':
    adaptive_boundary_selection()
    main()




