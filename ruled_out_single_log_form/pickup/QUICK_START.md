# Quick Start Guide

## Environment Setup

This project uses SageMath via Docker. The wrapper script `sage.ps1` handles execution.

```powershell
# From project root (C:\Users\zacha\physics)
.\sage.ps1 <script.sage>
```

## Verify Installation

```powershell
.\sage.ps1 tests/test_channel_identities.sage
# Expected: "SUCCESS: 50 samples verified for channel identities."
```

## Run All Tests

```powershell
.\sage.ps1 tools/run_all_tests.sage
```

## Key Entry Points

| Task | Command |
|------|---------|
| Test Hodges invariance | `.\sage.ps1 tests/test_hodges_invariance.sage` |
| Test factorization | `.\sage.ps1 tests/test_factorization_oracle.sage` |
| Run DCP geometry search | `.\sage.ps1 54.sage` |
| Check amplitude weights | `.\sage.ps1 tools/check_weights.sage` |

## File Organization

```
physics/
├── src/                    # Core library
│   ├── hodges.sage         # Hodges amplitude
│   ├── klt.sage            # KLT amplitude  
│   ├── spinor_sampling.sage# Kinematics generator
│   ├── oracle_gravity_mhv6.sage # Unified interface
│   └── boundary_sampler.sage   # BCFW shifts
├── tests/                  # Validation scripts
├── tools/                  # Utilities
├── results/                # Output reports
├── pickup/                 # Handoff documents (you are here)
└── 54.sage                 # Main DCP search engine
```

## Current Task

**Phase 5**: Prove uniqueness of gravity form via geometric constraints.

See `AGENT_HANDOFF.md` for full context.

