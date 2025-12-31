#!/usr/bin/env sage
# =============================================================================
# FAST GRAVITY POSITIVE GEOMETRY SEARCH - OPTIMIZED FOR 12 CORES
# =============================================================================
# Optimizations:
# - Uses all 12 CPU cores
# - Optimized memory usage
# - Fast PARI kernel
# - Strategic boundary selection
# - Early termination when dim=1 found
# =============================================================================

# Set configuration BEFORE loading 54.sage
FORCE_INVARIANT_MODE = 'S3xS3Z2'  # Start with S3xS3Z2 (dim=26, more promising than S6)

# Load main script
load('54.sage')

# =============================================================================
# OPTIMIZATION OVERRIDES FOR 12-CORE SYSTEM
# =============================================================================
import os

# Use all 12 cores (minus 1 for system)
_TOTAL_CORES = 12  # Explicitly set for your system
SCAN_WORKERS = 11  # Use 11 cores for parallel processing

# Optimize chunk size for better parallelization
MP_CHUNKSIZE = 4  # Larger chunks = less overhead

# Increase trials slightly for better chart quality
TOTAL_TRIALS = 75000  # Balance between speed and quality

# Use all 10 boundaries
INTERSECT_BOUNDARY_MODE = "ALL_3x3"

# Enable early termination
INTERSECT_TARGET_DIM = 1
INTERSECT_STOP_ON_EMPTY = False  # Don't stop if one boundary fails, try others

# Multi-strategy: try S3xS3Z2 first, then S3xS3, then S6
MULTI_STRATEGY_SEARCH = True
STRATEGY_INVARIANT_MODES = ['S3xS3Z2', 'S3xS3', 'S6']  # Start with most promising
STRATEGY_STOP_ON_DIM1 = True

# Optimize memory
CACHE_PRECOMPUTE = True
CACHE_CHARTS = True
gc.collect()  # Clean up before starting

# =============================================================================
# MAIN EXECUTION
# =============================================================================
log("\n" + "="*70)
log("FAST GRAVITY SEARCH - OPTIMIZED FOR 12 CORES")
log("="*70)
log(f"CPU Cores: {_TOTAL_CORES}")
log(f"Workers: {SCAN_WORKERS}")
log(f"Chunk Size: {MP_CHUNKSIZE}")
log(f"Trials: {TOTAL_TRIALS}")
log(f"Strategy: {STRATEGY_INVARIANT_MODES}")
log("="*70 + "\n")

# Run the search
main()


