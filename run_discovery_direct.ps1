# Direct discovery runner - runs phases sequentially
Write-Host "Starting Gravity Positive Geometry Discovery..."
Write-Host ""

$ErrorActionPreference = "Continue"

# Phase 1
Write-Host "="*70
Write-Host "PHASE 1: Pushforward Diagnostic"
Write-Host "="*70
docker run --rm -v "${PWD}:/workspace" sage-cursor bash -c "cd /workspace && sage src/discovery/phase1_pushforward.sage" 2>&1 | Tee-Object -Append discovery_results/phase1.log

# Phase 2  
Write-Host ""
Write-Host "="*70
Write-Host "PHASE 2: BCFW Amplituhedron"
Write-Host "="*70
docker run --rm -v "${PWD}:/workspace" sage-cursor bash -c "cd /workspace && sage src/discovery/phase2_bcfw.sage" 2>&1 | Tee-Object -Append discovery_results/phase2.log

# Phase 3
Write-Host ""
Write-Host "="*70
Write-Host "PHASE 3: Positivity Search"
Write-Host "="*70
docker run --rm -v "${PWD}:/workspace" sage-cursor bash -c "cd /workspace && sage src/discovery/phase3_positivity.sage" 2>&1 | Tee-Object -Append discovery_results/phase3.log

# Verdict
Write-Host ""
Write-Host "="*70
Write-Host "GENERATING VERDICT"
Write-Host "="*70
docker run --rm -v "${PWD}:/workspace" sage-cursor bash -c "cd /workspace && sage src/discovery/verdict.sage" 2>&1 | Tee-Object -Append discovery_results/verdict.log

Write-Host ""
Write-Host "="*70
Write-Host "DISCOVERY COMPLETE"
Write-Host "="*70
Write-Host "Check discovery_results/VERDICT.md for final conclusion"


