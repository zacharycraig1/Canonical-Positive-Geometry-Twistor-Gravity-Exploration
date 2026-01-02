# 🚀 START HERE: 6-Point MHV Gravity Positive Geometry Search

## Quick Start

**Run the search:**
```bash
sage 54.sage
```

Or for iterative search with auto-analysis:
```bash
sage run_iterative_search.sage
```

## What This Does

Searches for a **dimension-1 positive geometry** that describes 6-point MHV gravity amplitudes by:

1. Building the matroid structure for 6-point amplitudes
2. Computing invariant spaces under symmetry groups
3. Intersecting across boundaries to reduce dimension
4. Finding the unique (dim=1) candidate

## Current Configuration

- **Multi-strategy enabled**: Tries S3xS3 → S3xS3Z2 → S6
- **All 10 boundaries**: Uses all distinct 3|3 splits
- **Checkpointing**: Saves progress after each boundary
- **Target**: Dimension 1 (unique up to scale)

## Expected Runtime

- Per boundary: 10-60 seconds (chart search + scan)
- Total for 10 boundaries: 5-15 minutes
- Multi-strategy: Up to 3x longer if needed

## Results

Results saved to: `runs/YYYY-MM-DD_HHMMSS_seed42_modeS3xS3/`

Key files:
- `hit_report_intersection.json` - Full results
- `dim1_candidate_verification.json` - If dim=1 found
- `intersection_checkpoint.sobj` - Resume point

## If Interrupted

Set `RESUME_FROM_CHECKPOINT = True` and run again - it will resume from last boundary.

## Documentation

- `ITERATION_LOG.md` - Track results here
- `OPTIMIZATION_NOTES.md` - Performance details
- `PERFORMANCE_GUIDE.md` - Tuning guide
- `EXECUTION_READY.md` - Full setup details

---

**Ready to find the gravity positive geometry!** 🎯













