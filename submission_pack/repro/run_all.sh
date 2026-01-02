#!/bin/bash
set -e

echo "Running Full Verification Suite..."
echo "=================================="

echo ""
echo "[1/9] Reference Independence (F1.1)"
sage -python src/scripts/phaseF1_reference_independence.py

echo ""
echo "[2/9] Deletion Set Independence (F1.2)"
sage -python src/scripts/phaseF2_deletion_set_independence.py

echo ""
echo "[3/9] Pole Order Audit (F1.3)"
sage -python src/scripts/phaseF3_exact_pole_orders.py

echo ""
echo "[4/9] Forest Expansion (F2.1)"
sage -python src/scripts/phaseF4_all_minors_forest_expansion.py

echo ""
echo "[5/9] Newton Polytopes (F3.1)"
sage -python src/scripts/phaseF5_newton_polytopes.py

echo ""
echo "[6/9] n=7 Identity Verification (F4.1)"
sage -python src/scripts/phaseF6_n7_verification.py

echo ""
echo "[7/9] n=7 Valuations (F4.2)"
sage -python src/scripts/phaseF7_n7_valuations.py

echo ""
echo "[8/9] MTT Consistency Check"
sage tests/test_oracle_match.sage

echo ""
echo "[9/9] n=7 Sign Rule Verification"
sage src/signed_geometry/generalize_n7.sage

echo ""
echo "=================================="
echo "All verifications passed!"
