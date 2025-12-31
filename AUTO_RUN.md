# Auto-Run Configuration

## Optimizations Applied

1. **Reduced search parameters:**
   - TOTAL_TRIALS: 100000 → 50000 (faster iteration)
   - MAX_SOLUTIONS: 300 → 100 (faster processing)
   - VALIDATION_SHUFFLES: 2 → 1 (faster validation)

2. **Chart caching enabled:**
   - Caches charts per boundary
   - Reuses when same boundary/dimension encountered
   - Significant speedup on repeated runs

3. **Multi-strategy enabled:**
   - Automatically tries S3xS3 → S3xS3Z2 → S6
   - Stops when dim=1 found

4. **Auto-analysis:**
   - `analyze_and_optimize.py` analyzes results
   - Suggests next parameters automatically
   - `run_optimized_search.sage` applies suggestions

## To Run

**Best option (auto-iterative):**
```bash
sage run_optimized_search.sage
```

**Quick test:**
```bash
sage quick_search.sage
```

**Full search:**
```bash
sage 54.sage
```

## What Happens

1. Runs search with current parameters
2. Analyzes results automatically
3. Adjusts parameters based on results
4. Repeats until dim=1 found
5. Stops and reports success

## Progress Tracking

- Check `ITERATION_LOG.md` for results
- Each run saves to `runs/` directory
- Suggestions saved in each run directory









