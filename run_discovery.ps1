# Gravity Positive Geometry Discovery Runner
# This script runs the full discovery pipeline

Write-Host "="*70
Write-Host "GRAVITY POSITIVE GEOMETRY DISCOVERY"
Write-Host "="*70
Write-Host ""

# Check if Docker is available
try {
    docker --version | Out-Null
    Write-Host "[OK] Docker found"
} catch {
    Write-Host "[ERROR] Docker not found. Please install Docker."
    exit 1
}

# Check if Sage image exists
$imageExists = docker images sage-cursor -q
if (-not $imageExists) {
    Write-Host "[INFO] Building Docker image..."
    docker build -t sage-cursor .
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERROR] Docker build failed"
        exit 1
    }
}

# Create results directory
New-Item -ItemType Directory -Force -Path discovery_results | Out-Null

Write-Host "[INFO] Starting discovery pipeline..."
Write-Host "[INFO] This will take 8-10 hours. Results will be saved to discovery_results/"
Write-Host ""

# Run master orchestrator
$logFile = "discovery_results/discovery_$(Get-Date -Format 'yyyyMMdd_HHmmss').log"
Write-Host "[INFO] Logging to: $logFile"
Write-Host ""

.\sage.ps1 src/discovery/master_orchestrator.sage 2>&1 | Tee-Object $logFile

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "[SUCCESS] Discovery pipeline completed!"
    Write-Host "[INFO] Generating verdict..."
    
    # Generate verdict
    .\sage.ps1 src/discovery/verdict.sage 2>&1 | Tee-Object -Append $logFile
    
    Write-Host ""
    Write-Host "="*70
    Write-Host "DISCOVERY COMPLETE"
    Write-Host "="*70
    Write-Host "Results: discovery_results/"
    Write-Host "Verdict: discovery_results/VERDICT.md"
    Write-Host "="*70
} else {
    Write-Host ""
    Write-Host "[ERROR] Discovery pipeline failed. Check log: $logFile"
    exit 1
}


