#!/bin/bash
set -e

echo "Running Phase F Verification Suite..."

echo "[1/7] Reference Independence (F1.1)"
sage -python src/scripts/phaseF1_reference_independence.py

echo "[2/7] Deletion Set Independence (F1.2)"
sage -python src/scripts/phaseF2_deletion_set_independence.py

echo "[3/7] Pole Order Audit (F1.3)"
sage -python src/scripts/phaseF3_exact_pole_orders.py

echo "[4/7] Forest Expansion (F2.1)"
sage -python src/scripts/phaseF4_all_minors_forest_expansion.py

echo "[5/7] Newton Polytopes (F3.1)"
sage -python src/scripts/phaseF5_newton_polytopes.py

echo "[6/7] n=7 Verification (F4.1)"
sage -python src/scripts/phaseF6_n7_verification.py

echo "[7/7] n=7 Valuations (F4.2)"
sage -python src/scripts/phaseF7_n7_valuations.py

echo "All verifications passed!"







