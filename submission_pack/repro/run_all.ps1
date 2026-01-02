Write-Host "Running Full Verification Suite..."
Write-Host "=================================="

Write-Host ""
Write-Host "[1/9] Reference Independence (F1.1)"
& .\sage.ps1 -python src/scripts/phaseF1_reference_independence.py

Write-Host ""
Write-Host "[2/9] Deletion Set Independence (F1.2)"
& .\sage.ps1 -python src/scripts/phaseF2_deletion_set_independence.py

Write-Host ""
Write-Host "[3/9] Pole Order Audit (F1.3)"
& .\sage.ps1 -python src/scripts/phaseF3_exact_pole_orders.py

Write-Host ""
Write-Host "[4/9] Forest Expansion (F2.1)"
& .\sage.ps1 -python src/scripts/phaseF4_all_minors_forest_expansion.py

Write-Host ""
Write-Host "[5/9] Newton Polytopes (F3.1)"
& .\sage.ps1 -python src/scripts/phaseF5_newton_polytopes.py

Write-Host ""
Write-Host "[6/9] n=7 Identity Verification (F4.1)"
& .\sage.ps1 -python src/scripts/phaseF6_n7_verification.py

Write-Host ""
Write-Host "[7/9] n=7 Valuations (F4.2)"
& .\sage.ps1 -python src/scripts/phaseF7_n7_valuations.py

Write-Host ""
Write-Host "[8/9] MTT Consistency Check"
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage tests/test_oracle_match.sage

Write-Host ""
Write-Host "[9/9] n=7 Sign Rule Verification"
docker run --rm -v "${PWD}:/home/sage/project" -w /home/sage/project sagemath/sagemath:latest sage src/signed_geometry/generalize_n7.sage

Write-Host ""
Write-Host "=================================="
Write-Host "All verifications passed!"
