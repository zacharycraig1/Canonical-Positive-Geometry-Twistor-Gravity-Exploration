# Performance Guide

## Quick Performance Tips

### If Running Too Slowly:

1. **Enable Checkpointing** (already enabled by default)
   - Progress is saved after each boundary
   - Can resume if interrupted
   - No performance cost, only saves time on resume

2. **Reduce Search Trials** (if accuracy allows)
   ```python
   TOTAL_TRIALS = 50000  # Reduce from 100000
   ```
   - May find slightly smaller charts
   - Still finds valid charts, just not always optimal

3. **Use Fewer Boundaries** (if dim=1 not needed immediately)
   ```python
   INTERSECT_BOUNDARY_MODE = "CUSTOM"
   INTERSECT_BOUNDARIES = [(1,2,3), (1,2,4), (1,3,4)]  # Start with 3
   ```
   - Faster, but may not reach dim=1
   - Can add more boundaries later

4. **Single Strategy Mode** (faster than multi-strategy)
   ```python
   MULTI_STRATEGY_SEARCH = False  # Already default
   ```

5. **Check Worker Count**
   ```python
   SCAN_WORKERS = max(1, os.cpu_count() - 1)  # Auto-scales
   ```
   - Should match your CPU cores
   - More workers = faster parallel operations

## Performance Monitoring

The script reports:
- `Chart search: X charts found in Y.Ys` - Time per boundary
- `Intersection complete in Y.Ys` - Total intersection time
- `rate=X/s` - Processing rates

**Watch for:**
- Chart search > 60s → Consider reducing `TOTAL_TRIALS`
- Nullspace computation slow → Already optimized with PARI
- Memory issues → Checkpointing helps

## Resuming from Checkpoint

If interrupted:
```python
RESUME_FROM_CHECKPOINT = True
# Then run normally
sage 54.sage
```

The script will:
1. Detect checkpoint file
2. Resume from last boundary
3. Continue from there

## Accuracy vs Speed Trade-offs

### Safe Optimizations (No Accuracy Loss):
- ✅ Checkpointing
- ✅ Early termination (when enough solutions found)
- ✅ PARI kernel (mathematically equivalent)
- ✅ Parallel processing
- ✅ Modular prefilter

### Potential Trade-offs (Use with Caution):
- ⚠️ Reducing `TOTAL_TRIALS` - May miss optimal charts
- ⚠️ Fewer boundaries - May not reach dim=1
- ⚠️ Single strategy - May miss better invariant mode

## Expected Performance

**Typical run times:**
- Matroid/OS3 construction: < 1s
- Invariant computation: 1-5s
- Per boundary:
  - Chart search: 10-60s (depends on TOTAL_TRIALS)
  - Chart scan: 1-30s (depends on chart size)
- Total for 10 boundaries: 5-15 minutes

**With optimizations:**
- Early termination: Can cut time in half if good charts found early
- Checkpointing: Prevents rework on interruption
- Parallel workers: 2-4x speedup on multi-core systems

## Troubleshooting Slow Runs

1. **Check CPU usage**
   - Should be high during parallel operations
   - If low, may be I/O bound (check disk speed)

2. **Check memory**
   - Large matrices can use significant RAM
   - Checkpointing helps by allowing restarts

3. **Monitor progress**
   - Watch boundary-by-boundary progress
   - Identify which boundary is slow

4. **Adjust parameters**
   - Reduce `TOTAL_TRIALS` if chart search is slow
   - Reduce boundaries if intersection is slow
   - Increase workers if CPU underutilized













