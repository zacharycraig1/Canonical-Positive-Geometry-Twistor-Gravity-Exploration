# Running the Dimension-1 Candidate Search

## Quick Start

```bash
sage 54.sage
```

## Configuration Options

### Basic Settings (in `54.sage`)

**Single Strategy Mode (default):**
```python
MULTI_STRATEGY_SEARCH = False
INVARIANT_MODE = 'S3xS3'  # Options: 'S6', 'S3xS3', 'S3xS3Z2'
```

**Multi-Strategy Mode (recommended for exploration):**
```python
MULTI_STRATEGY_SEARCH = True
STRATEGY_INVARIANT_MODES = ['S3xS3', 'S3xS3Z2', 'S6']
STRATEGY_STOP_ON_DIM1 = True  # Stop when dim=1 found
```

### Search Parameters

```python
TOTAL_TRIALS = 100000        # Chart search trials per seed
NUM_SEEDS = 1                # Number of random seeds
BASE_SEED = 42               # Starting seed
INTERSECT_TARGET_DIM = 1     # Target dimension
INTERSECT_BOUNDARY_MODE = "ALL_3x3"  # Use all 10 boundaries
```

### Performance Tuning

```python
SCAN_WORKERS = max(1, os.cpu_count() - 1)  # Auto-scales to cores
USE_PARI_KERNEL = True                      # 10-100x speedup
MODP_PREFILTER = True                       # Modular prefilter
```

## What the Script Does

1. **Builds the matroid M6** for 6-point amplitudes
2. **Constructs OS3 space** (outer space)
3. **Computes invariants** under symmetry group (S3×S3, S6, etc.)
4. **Intersects across boundaries** to reduce dimension
5. **Stops when dim ≤ 1** or intersection becomes empty
6. **Verifies dim=1 candidates** automatically
7. **Saves all artifacts** to `runs/` directory

## Output

### When dim=1 Candidate Found:
```
[CANDIDATE] Candidate direction found (unique up to scale on this run).
[SUCCESS] Dimension-1 candidate verified and ready for Hodge structure check!
```

### Artifacts Saved:
- `runs/YYYY-MM-DD_HHMMSS_seed42_modeS3xS3/`
  - `config.json` - Configuration used
  - `hit_report_intersection.json` - Full results
  - `dim1_candidate_verification.json` - Verification results (if dim=1)
  - `candidate_space_basis_*.sobj` - Basis vectors
  - `alpha0_first_basis_vector.sobj` - The candidate

## Troubleshooting

### If dim > 1:
- Try different `INVARIANT_MODE` (S6, S3xS3, S3xS3Z2)
- Enable `MULTI_STRATEGY_SEARCH = True` to try all modes
- Increase `TOTAL_TRIALS` for better chart search
- Try different `BASE_SEED` values

### If intersection becomes empty:
- Try fewer boundaries (use `INTERSECT_BOUNDARY_MODE = "CUSTOM"`)
- Relax constraints
- Try different invariant modes

### Performance Issues:
- Ensure `USE_PARI_KERNEL = True` (default)
- Check `SCAN_WORKERS` matches your CPU cores
- Enable `CACHE_PRECOMPUTE = True` (default)

## Next Steps After Finding dim=1

1. Review `dim1_candidate_verification.json`
2. Check `alpha0_first_basis_vector.sobj` for the candidate
3. Implement Hodge structure computation
4. Verify against expected values for 6pt MHV gravity












