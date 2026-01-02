# Gravity Positive Geometry Discovery - Robust Runner
# This script runs the discovery pipeline with proper error handling

$ErrorActionPreference = "Continue"
$ProgressPreference = "SilentlyContinue"

Write-Host "="*70
Write-Host "GRAVITY POSITIVE GEOMETRY DISCOVERY"
Write-Host "="*70
Write-Host ""

# Ensure results directory exists
New-Item -ItemType Directory -Force -Path discovery_results | Out-Null

# Function to run a phase
function Run-Phase {
    param($PhaseName, $ScriptPath)
    
    Write-Host ""
    Write-Host "="*70
    Write-Host "PHASE: $PhaseName"
    Write-Host "="*70
    Write-Host "Running: $ScriptPath"
    Write-Host ""
    
    $logFile = "discovery_results/${PhaseName}_$(Get-Date -Format 'yyyyMMdd_HHmmss').log"
    
    $result = docker run --rm `
        -v "${PWD}:/workspace" `
        sage-cursor `
        bash -c "cd /workspace && sage $ScriptPath" `
        2>&1 | Tee-Object -FilePath $logFile
    
    $exitCode = $LASTEXITCODE
    
    if ($exitCode -eq 0) {
        Write-Host "[SUCCESS] $PhaseName completed"
    } else {
        Write-Host "[WARNING] $PhaseName exited with code $exitCode (check log: $logFile)"
    }
    
    return $exitCode
}

# Run Phase 1
$phase1Exit = Run-Phase "Phase1_Pushforward" "src/discovery/phase1_pushforward.sage"

# Run Phase 2  
$phase2Exit = Run-Phase "Phase2_BCFW" "src/discovery/phase2_bcfw.sage"

# Run Phase 3
$phase3Exit = Run-Phase "Phase3_Positivity" "src/discovery/phase3_positivity.sage"

# Generate Verdict
Write-Host ""
Write-Host "="*70
Write-Host "GENERATING VERDICT"
Write-Host "="*70

$verdictExit = Run-Phase "Verdict" "src/discovery/verdict.sage"

# Final summary
Write-Host ""
Write-Host "="*70
Write-Host "DISCOVERY PIPELINE COMPLETE"
Write-Host "="*70
Write-Host "Phase 1: $(if ($phase1Exit -eq 0) { 'SUCCESS' } else { 'FAILED' })"
Write-Host "Phase 2: $(if ($phase2Exit -eq 0) { 'SUCCESS' } else { 'FAILED' })"
Write-Host "Phase 3: $(if ($phase3Exit -eq 0) { 'SUCCESS' } else { 'FAILED' })"
Write-Host "Verdict:  $(if ($verdictExit -eq 0) { 'SUCCESS' } else { 'FAILED' })"
Write-Host ""
Write-Host "Results: discovery_results/"
Write-Host "Verdict: discovery_results/VERDICT.md"
Write-Host "="*70


