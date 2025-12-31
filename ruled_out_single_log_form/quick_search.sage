#!/usr/bin/env sage
# Quick optimized search - reduced parameters for faster iteration
# Run with: sage quick_search.sage

load('54.sage')

# Override for quick search
TOTAL_TRIALS = 30000  # Even faster
MAX_SOLUTIONS = 50
VALIDATION_SHUFFLES = 0  # Skip validation for speed
VALIDATION_SEED_OFFSETS = []
RUN_CANDIDATE_VALIDATION = False  # Skip for speed

# But keep multi-strategy
MULTI_STRATEGY_SEARCH = True
STRATEGY_STOP_ON_DIM1 = True

if __name__ == '__main__':
    main()




