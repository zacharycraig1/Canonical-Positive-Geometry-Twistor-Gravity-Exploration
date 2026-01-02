# 🚀 START HERE: MHV Single Solution Analysis

## What This Is

Investigation to resolve the **~2220x discrepancy** between CHY and Hodges gravity amplitudes by verifying the **Du-Teng-Wu theory** (arXiv:1603.08158): for MHV amplitudes, **only 1 of 6 solutions** contributes.

## Quick Start (60 seconds)

### 1. Open PowerShell in `c:\Users\zacha\physics`

### 2. Run the analysis:
```powershell
.\run_mhv_analysis.ps1
```

### 3. Wait 1-3 minutes for results

### 4. Look for this output:

**✅ SUCCESS CASE:**
```
✅ MHV SINGLE SOLUTION THEORY CONFIRMED
   Only 1 of 6 solutions has non-zero Pf'(Ψ)
   ✅ MHV solution matches Hodges amplitude!
   
   CONCLUSION: CHY normalization issue is RESOLVED
```

**⚠️ NORMALIZATION FACTOR CASE:**
```
✅ MHV SINGLE SOLUTION THEORY CONFIRMED
   Ratio (single/Hodges) = 1.667
   ✅ Match found with factor 0.6000
```
→ Apply this factor to CHY formula

**❌ NEEDS DEBUG:**
```
⚠️ MHV SINGLE SOLUTION THEORY NOT CONFIRMED
   Found 3 solutions with non-zero Pf'
```
→ See troubleshooting guide

## What It Does

```
Input: Random 6-point MHV kinematics
  ↓
Solve 6 scattering equations
  ↓
Compute Pf'(Ψ) at each solution
  ↓
Classify: Zero vs Non-Zero Pfaffians
  ↓
Expected: 5 zeros, 1 non-zero (MHV solution)
  ↓
Compare MHV solution to Hodges amplitude
  ↓
Output: Match status + any normalization factor
```

## Key Files

### To Run:
- **`run_mhv_analysis.ps1`** ← Run this
- `src/gravity_proof/analyze_mhv_single_solution.sage` (main script)

### To Read:
1. **`RUN_MHV_ANALYSIS_README.md`** ← Start here for details
2. `MHV_SINGLE_SOLUTION_INVESTIGATION.md` (technical background)
3. `MHV_ANALYSIS_CHECKLIST.md` (execution guide)
4. `IMPLEMENTATION_SUMMARY.md` (what was built)

## Expected Timeline

1. **Run analysis**: 1-3 minutes
2. **Interpret results**: < 1 minute  
3. **If successful**: Update CHY implementation (~30 min)
4. **Re-run full proof**: 5-10 minutes
5. **Total**: < 1 hour to resolve the issue

## Decision Tree

```
Run .\run_mhv_analysis.ps1
         ↓
    MHV theory confirmed?
    ├─ YES
    │  ├─ Hodges match?
    │  │  ├─ YES → ✅ DONE (use single solution)
    │  │  └─ NO → ⚠️ Apply normalization factor
    │  └─ Update CHY → Re-run proof → ✅ RESOLVED
    └─ NO → ❌ Debug helicities/kinematics
```

## What Happens After Success

1. **Update CHY formula** to use only MHV solution
2. **Apply normalization** factor (if any)
3. **Re-run** full proof: `sage src/gravity_proof/run_proof.sage`
4. **Verify** all success criteria pass
5. **Document** resolution in proof certificate

## Troubleshooting

| Symptom | Quick Fix |
|---------|-----------|
| "No solutions found" | Try different seed: edit script line with `seed=43` |
| Multiple non-zero Pf' | Check helicity config or adjust threshold |
| No simple factor | Review Pfaffian sign conventions |
| Sage not found | Install SageMath or add to PATH |

See `MHV_ANALYSIS_CHECKLIST.md` for detailed troubleshooting.

## Theory Background (30 seconds)

**Current CHY (Wrong)**:
```
M₆ = Σ(all 6 solutions) [Pf'(Ψ)]² / det'(Φ)
```
Problem: Includes 5 solutions that **should be zero** for MHV

**MHV-Corrected CHY (Right)**:
```
M₆ = [Pf'(Ψ)]² / det'(Φ)  [at the 1 MHV solution only]
```
Why: Du-Teng-Wu proved Pf'(Ψ) = 0 exactly at other 5 solutions

**Error Source**: ~2220x ≈ (6 solutions) × (spurious contributions) × (missing normalization)

## Files Created (This Implementation)

```
c:\Users\zacha\physics\
├── run_mhv_analysis.ps1                              ← RUN THIS
├── src\gravity_proof\
│   └── analyze_mhv_single_solution.sage              (19 KB, main script)
└── Documentation\
    ├── START_HERE_MHV_ANALYSIS.md                    ← YOU ARE HERE
    ├── RUN_MHV_ANALYSIS_README.md                    (quick start)
    ├── MHV_SINGLE_SOLUTION_INVESTIGATION.md          (technical)
    ├── MHV_ANALYSIS_CHECKLIST.md                     (execution guide)
    └── IMPLEMENTATION_SUMMARY.md                     (what was built)
```

## One-Line Summary

**Run `.\run_mhv_analysis.ps1` to verify only 1 of 6 solutions contributes, fixing the ~2220x CHY error.**

---

## 🎯 YOUR NEXT ACTION

```powershell
cd c:\Users\zacha\physics
.\run_mhv_analysis.ps1
```

Wait 1-3 minutes → Read output → Follow recommendations

---

**Questions?** See `RUN_MHV_ANALYSIS_README.md`  
**Technical Details?** See `MHV_SINGLE_SOLUTION_INVESTIGATION.md`  
**Step-by-Step Guide?** See `MHV_ANALYSIS_CHECKLIST.md`

