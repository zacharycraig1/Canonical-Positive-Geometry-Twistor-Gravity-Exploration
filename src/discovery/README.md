# Gravity Positive Geometry Discovery Pipeline

This directory contains the automated discovery pipeline to definitively determine whether MHV gravity amplitudes can be expressed as the canonical form of a positive geometry.

## Quick Start

```bash
# Run the full discovery pipeline (8-10 hours)
.\sage.ps1 src/discovery/master_orchestrator.sage 2>&1 | Tee-Object discovery.log

# Generate verdict from results
.\sage.ps1 src/discovery/verdict.sage
```

## Pipeline Overview

The discovery runs three parallel phases:

1. **Phase 1: Pushforward Diagnostic** (1 hour, 4 cores)
   - Tests 1000 kinematic points
   - Determines if pushforward/Hodges ratio is constant
   - If constant → normalization fix possible

2. **Phase 2: BCFW Amplituhedron** (4 hours, 8 cores)
   - Tests 10,000 kinematic points
   - Tests multiple BCFW formula variants
   - Finds formula that matches Hodges

3. **Phase 3: Positivity Search** (3 hours, 4 cores)
   - Brute force search: 10,000 random points
   - Optimization search: 100 restarts
   - Finds or proves non-existence of positive region

## Files

- `master_orchestrator.sage` - Main entry point, runs all phases
- `phase1_pushforward.sage` - Pushforward diagnostic
- `phase2_bcfw.sage` - BCFW formula testing
- `phase3_positivity.sage` - Positivity region search
- `verdict.sage` - Evidence collection and verdict generation
- `utils.sage` - Shared utilities

## Output

Results are saved to `discovery_results/`:

- `phase1_result.json` - Phase 1 results
- `phase2_result.json` - Phase 2 results
- `phase3_result.json` - Phase 3 results
- `all_phases_results.json` - Combined results
- `final_evidence.json` - Collected evidence
- `VERDICT.md` - Final verdict report
- `discovery.log` - Execution log

## Checkpointing

All phases support checkpointing - if interrupted, they resume from the last checkpoint.

## Running Individual Phases

```bash
# Phase 1 only
.\sage.ps1 src/discovery/phase1_pushforward.sage

# Phase 2 only
.\sage.ps1 src/discovery/phase2_bcfw.sage

# Phase 3 only
.\sage.ps1 src/discovery/phase3_positivity.sage
```

## Expected Runtime

- Full pipeline: 8-10 hours
- Phase 1: ~1 hour
- Phase 2: ~4 hours
- Phase 3: ~3 hours
- Verdict: ~30 minutes

## System Requirements

- 16GB RAM (recommended)
- 16 CPU cores (will use up to 16)
- Docker with SageMath image
- Python 3 with multiprocessing support


