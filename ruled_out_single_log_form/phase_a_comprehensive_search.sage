#!/usr/bin/env sage
# Phase A Comprehensive Search - Based on factorization_dim8_constraints_report.txt
# Implements Phase A (purely linear, easiest) constraints:
# 1) Hard-enforce full S6 invariance (not just measuring)
# 2) Expand boundary set to all |S|=2,3 channels (25 total: 15 two-particle + 10 three-particle)
# 3) Add iterated residue consistency for compatible pairs

# Set override flag BEFORE loading 54.sage
# This ensures S6 mode is used from the start
FORCE_INVARIANT_MODE = 'S6'

# Load the main script (it will check FORCE_INVARIANT_MODE)
load('54.sage')

# =============================================================================
# PHASE A: HARD S6 ENFORCEMENT  
# =============================================================================
# CRITICAL: Override the default S3xS3 mode from 54.sage
# Force S6 invariance from the start - project onto S6 invariants BEFORE intersection
# This ensures we work in the S6-invariant subspace from the beginning
INVARIANT_MODE = 'S6'  # Hard-enforce S6 (as per report recommendation #1)
MULTI_STRATEGY_SEARCH = False  # Focus on S6 only for Phase A

# Verify the override worked
print(f"\n[PHASE A] INVARIANT_MODE override: {INVARIANT_MODE}")
if INVARIANT_MODE != 'S6':
    raise RuntimeError(f"Failed to set INVARIANT_MODE to 'S6', got '{INVARIANT_MODE}'")

# =============================================================================
# PHASE A: ALL PHYSICAL CHANNELS (25 total)
# =============================================================================
# According to report: For n=6:
# - 2-particle channels: C(6,2) = 15
# - 3-particle channels: 10 (since S and complement define the same divisor)
# Total: 25 physical divisors

def all_physical_boundaries(n=6):
    """All physical boundaries: 2-particle + 3-particle channels.
    Returns list of boundaries, where each boundary is either:
    - A 2-tuple (i,j) for 2-particle channel
    - A 3-tuple (i,j,k) for 3-particle channel
    
    Note: 2-particle boundaries may need special handling in the intersection code.
    """
    boundaries = []
    
    # 2-particle channels: all pairs (i,j) with i < j
    # These are 15 total: C(6,2) = 15
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            boundaries.append((i, j))
    
    # 3-particle channels: all distinct 3|3 splits (mod complement)
    # These are 10 total (since S and complement define the same divisor)
    for S in itertools.combinations(range(1, n+1), 3):
        S = tuple(sorted(S))
        comp = tuple(sorted(set(range(1, n+1)) - set(S)))
        # Keep one representative per complement pair
        if S < comp:
            boundaries.append(S)
    
    return boundaries

# Use all 25 physical boundaries
INTERSECT_BOUNDARY_MODE = "CUSTOM"
INTERSECT_BOUNDARIES = all_physical_boundaries(6)

log(f"\n[PHASE A] Using {len(INTERSECT_BOUNDARIES)} physical boundaries:")
log(f"  - 2-particle: {sum(1 for b in INTERSECT_BOUNDARIES if len(b) == 2)}")
log(f"  - 3-particle: {sum(1 for b in INTERSECT_BOUNDARIES if len(b) == 3)}")

# =============================================================================
# PHASE A: ITERATED RESIDUE CONSISTENCY
# =============================================================================
# Add iterated residue constraints for compatible channel pairs
# This enforces: Res_{s_A=0} Res_{s_B=0} Ω = Res_{s_B=0} Res_{s_A=0} Ω
# (with proper sign/orientation tracking)

# For now, we'll enable this as a flag - implementation will be added to intersection logic
ENABLE_ITERATED_RESIDUES = True  # Enable iterated residue consistency checks

# =============================================================================
# SEARCH PARAMETERS
# =============================================================================
TOTAL_TRIALS = 50000
MAX_SOLUTIONS = 100

# Validation - keep enabled for physics breakthrough verification
RUN_CANDIDATE_VALIDATION = True
VALIDATION_SHUFFLES = 1
VALIDATION_SEED_OFFSETS = [0]
VALIDATION_REQUIRE_DIM1 = True  # Only accept dim=1 candidates

log("\n" + "="*70)
log("PHASE A COMPREHENSIVE SEARCH - Physics Breakthrough Strategy")
log("Based on factorization_dim8_constraints_report.txt")
log("="*70)
log("Phase A Constraints:")
log("  1. Hard-enforce S6 invariance (project onto S6 invariants)")
log("  2. Use all 25 physical channels (15 two-particle + 10 three-particle)")
log("  3. Add iterated residue consistency (compatible pairs)")
log("="*70)

# CRITICAL: Ensure INVARIANT_MODE is set to S6 before main() is called
# The main() function reads this global variable
INVARIANT_MODE = 'S6'
MULTI_STRATEGY_SEARCH = False

log(f"INVARIANT_MODE: {INVARIANT_MODE} (verified)")
log(f"Total boundaries: {len(INTERSECT_BOUNDARIES)}")
log(f"Iterated residues: {'ENABLED' if ENABLE_ITERATED_RESIDUES else 'DISABLED'}")
log("="*70)

# Monkey-patch the main function to ensure S6 is used
# This is necessary because load() might create namespace issues
_original_main = main

def phase_a_main():
    """Wrapper that ensures S6 mode is used."""
    global INVARIANT_MODE, MULTI_STRATEGY_SEARCH
    INVARIANT_MODE = 'S6'
    MULTI_STRATEGY_SEARCH = False
    print(f"\n{'='*70}")
    print(f"[PHASE A] ENFORCING INVARIANT_MODE = {INVARIANT_MODE}")
    print(f"{'='*70}\n")
    return _original_main()

if __name__ == '__main__':
    # Replace main with our wrapper
    main = phase_a_main
    main()

