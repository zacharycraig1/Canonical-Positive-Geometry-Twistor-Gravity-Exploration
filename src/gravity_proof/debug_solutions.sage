#!/usr/bin/env sage
"""
Debug script to examine scattering equation solutions in detail.
"""
from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

load("src/gravity_proof/scattering_solver.sage")
load("src/gravity_proof/psi_matrix.sage")

from src.kinematics.spinors import SpinorKinematics

print("="*80)
print("DEBUG: Examining Scattering Equation Solutions")
print("="*80)

# Use same seed as main analysis
kin = SpinorKinematics.random_rational(6, seed=42)
solver = ScatteringEquationSolver(kin)

print("\nSolving scattering equations...")
solutions = solver.solve_numerical()

print(f"\nFound {len(solutions)} solutions")
print("\n" + "="*80)
print("SOLUTION DETAILS")
print("="*80)

for i, sol in enumerate(solutions):
    z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
    
    # Real parts - handle both Sage types and Python floats
    def to_real_float(x):
        if hasattr(x, 'real') and callable(x.real):
            return float(x.real())
        elif hasattr(x, 'real'):
            return float(x.real)
        else:
            return float(x)
    
    z4_r = to_real_float(z4)
    z5_r = to_real_float(z5)
    z6_r = to_real_float(z6)
    
    print(f"\nSolution {i+1}:")
    print(f"  z4 = {z4_r:+.6f}")
    print(f"  z5 = {z5_r:+.6f}")
    print(f"  z6 = {z6_r:+.6f}")
    
    # Check various orderings
    orderings = [
        ("z4 > z5 > z6 > 0", z4_r > z5_r > z6_r > 0),
        ("z4 > z5 > z6", z4_r > z5_r > z6_r),
        ("z4 > z5, z4 > z6", z4_r > z5_r and z4_r > z6_r),
        ("All positive", z4_r > 0 and z5_r > 0 and z6_r > 0),
        ("z6 < z5 < z4", z6_r < z5_r < z4_r),
        ("0 < z6 < z5 < z4 < 1", 0 < z6_r < z5_r < z4_r < 1),
        ("Between 0 and 1", 0 < z4_r < 1 and 0 < z5_r < 1 and 0 < z6_r < 1),
    ]
    
    for desc, check in orderings:
        if check:
            print(f"  ✓ {desc}")
        else:
            print(f"  ✗ {desc}")
    
    # Compute Pfaffian at this point
    try:
        psi = PsiMatrixMHV(kin, z4_r, z5_r, z6_r)
        pf = psi.reduced_pfaffian_standard((0, 1))
        if pf == 0:
            pf = psi.reduced_pfaffian_delete_3_6()
        pf_val = float(pf) if hasattr(pf, '__float__') else complex(pf).real
        print(f"  Pf'(Psi) = {pf_val:.6e} ({'> 0' if pf_val > 0 else '< 0'})")
    except Exception as e:
        print(f"  Pf'(Psi) = ERROR: {e}")

# Summary of what orderings solutions satisfy
print("\n" + "="*80)
print("ORDERING ANALYSIS")
print("="*80)

# Test all 6 possible orderings of z4, z5, z6
orderings_count = {}
for i, sol in enumerate(solutions):
    z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
    z4_r = to_real_float(z4)
    z5_r = to_real_float(z5)
    z6_r = to_real_float(z6)
    
    # Determine ordering
    vals = [('z4', z4_r), ('z5', z5_r), ('z6', z6_r)]
    sorted_vals = sorted(vals, key=lambda x: x[1], reverse=True)
    ordering = ' > '.join([v[0] for v in sorted_vals])
    
    if ordering not in orderings_count:
        orderings_count[ordering] = []
    orderings_count[ordering].append(i+1)

print("\nSolutions by ordering:")
for ordering, sols in orderings_count.items():
    print(f"  {ordering}: solutions {sols}")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)

# Check if there's one ordering with all 6 solutions
if any(len(sols) == 6 for sols in orderings_count.values()):
    for ordering, sols in orderings_count.items():
        if len(sols) == 6:
            print(f"\n✓ All 6 solutions have ordering: {ordering}")
            print("  This should define the correct region R6!")
else:
    print("\n✗ Solutions have DIFFERENT orderings - no single ordering contains all 6")
    print("  This suggests the positive geometry is more complex than a simple ordering")
    print("\nPossible implications:")
    print("  1. The region R6 might be a union of multiple ordered regions")
    print("  2. Ordering constraints alone are insufficient")
    print("  3. Need additional positivity conditions beyond B1-B4")

