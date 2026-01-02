# =============================================================================
# CHY Sanity Checks Runner
# =============================================================================
# Runs all sanity checks to verify CHY implementation

Write-Host "============================================================"
Write-Host "CHY Implementation Sanity Checks"
Write-Host "============================================================"
Write-Host ""

# Ensure RESULTS directory exists
if (!(Test-Path "RESULTS")) {
    New-Item -ItemType Directory -Path "RESULTS" | Out-Null
}

# Test 1: n=4 CHY Test
Write-Host "============================================================"
Write-Host "TEST 1: n=4 CHY (Single Solution)"
Write-Host "============================================================"
Write-Host ""

docker run --rm -v ${PWD}:/workspace sage-cursor sage src/gravity_proof/test_n4_chy.sage 2>&1 | Tee-Object -FilePath "RESULTS/n4_chy_test.log"

Write-Host ""
Write-Host "============================================================"
Write-Host "TEST 2: BCFW vs Hodges Verification"
Write-Host "============================================================"
Write-Host ""

# Check if file exists
if (Test-Path "src/amplituhedron/verify_bcfw_hodges.sage") {
    docker run --rm -v ${PWD}:/workspace sage-cursor sage src/amplituhedron/verify_bcfw_hodges.sage 2>&1 | Tee-Object -FilePath "RESULTS/bcfw_hodges_test.log"
} else {
    Write-Host "SKIP: verify_bcfw_hodges.sage not found"
}

Write-Host ""
Write-Host "============================================================"
Write-Host "Sanity Checks Complete"
Write-Host "============================================================"
Write-Host "Logs saved to RESULTS/ directory"


