Write-Host "Running Phase F Verification Suite..."

Write-Host "[1/7] Reference Independence (F1.1)"
& .\sage.ps1 -python src/scripts/phaseF1_reference_independence.py

Write-Host "[2/7] Deletion Set Independence (F1.2)"
& .\sage.ps1 -python src/scripts/phaseF2_deletion_set_independence.py

Write-Host "[3/7] Pole Order Audit (F1.3)"
& .\sage.ps1 -python src/scripts/phaseF3_exact_pole_orders.py

Write-Host "[4/7] Forest Expansion (F2.1)"
& .\sage.ps1 -python src/scripts/phaseF4_all_minors_forest_expansion.py

Write-Host "[5/7] Newton Polytopes (F3.1)"
& .\sage.ps1 -python src/scripts/phaseF5_newton_polytopes.py

Write-Host "[6/7] n=7 Verification (F4.1)"
& .\sage.ps1 -python src/scripts/phaseF6_n7_verification.py

Write-Host "[7/7] n=7 Valuations (F4.2)"
& .\sage.ps1 -python src/scripts/phaseF7_n7_valuations.py

Write-Host "All verifications passed!"



