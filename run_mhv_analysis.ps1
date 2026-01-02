# PowerShell script to run MHV single-solution analysis
# Usage: .\run_mhv_analysis.ps1

Write-Host "================================================================" -ForegroundColor Cyan
Write-Host "MHV SINGLE SOLUTION ANALYSIS" -ForegroundColor Cyan  
Write-Host "Investigating CHY Normalization (Du-Teng-Wu arXiv:1603.08158)" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan

# Check if sage is available
try {
    sage --version | Out-Null
    Write-Host "`n✓ SageMath found" -ForegroundColor Green
} catch {
    Write-Host "`n✗ SageMath not found. Please install SageMath or add it to PATH." -ForegroundColor Red
    exit 1
}

# Run the analysis
Write-Host "`nRunning MHV single-solution analysis..." -ForegroundColor Yellow
Write-Host "(This may take several minutes)`n" -ForegroundColor Yellow

sage src\gravity_proof\analyze_mhv_single_solution.sage

$exitCode = $LASTEXITCODE

if ($exitCode -eq 0) {
    Write-Host "`n================================================================" -ForegroundColor Green
    Write-Host "ANALYSIS COMPLETED SUCCESSFULLY" -ForegroundColor Green
    Write-Host "================================================================" -ForegroundColor Green
} else {
    Write-Host "`n================================================================" -ForegroundColor Red
    Write-Host "ANALYSIS FAILED (exit code: $exitCode)" -ForegroundColor Red
    Write-Host "================================================================" -ForegroundColor Red
}

exit $exitCode

