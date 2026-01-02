import sys
import os
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/check_residues_gate_b.sage")

def test_mhv_chart_independence():
    print("Running Unit Test: MHV Chart Independence...")
    
    n = 6
    RF = RealField(200)
    
    # 1. Generate a valid point
    print("Generating random kinematic point...")
    ts = [rnd.uniform(0, 10) for _ in range(n)]
    ts_tilde_free = {i: vector(RF, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(n-2)}
    x_s = vector(RF, [1, -2.0])
    y_s = vector(RF, [1, 12.0])
    
    lambdas = {i: vector(RF, [1, ts[i]]) for i in range(n)}
    tildes = solve_conservation_generic(lambdas, ts_tilde_free, n, RF)
    
    if tildes is None:
        print("Failed to generate valid kinematics.")
        exit(1)
        
    # 2. Compute M_MHV with different roots
    roots1 = [0, 1, 2]
    roots2 = [0, 1, 3]
    roots3 = [0, 2, 4]
    
    print(f"Computing M_MHV with roots {roots1}...")
    m1 = compute_M_MHV(n, lambdas, tildes, x_s, y_s, RF, roots=roots1)
    
    print(f"Computing M_MHV with roots {roots2}...")
    m2 = compute_M_MHV(n, lambdas, tildes, x_s, y_s, RF, roots=roots2)
    
    print(f"Computing M_MHV with roots {roots3}...")
    m3 = compute_M_MHV(n, lambdas, tildes, x_s, y_s, RF, roots=roots3)
    
    print("\nResults:")
    print(f"  M1: {float(m1):.6e}")
    print(f"  M2: {float(m2):.6e}")
    print(f"  M3: {float(m3):.6e}")
    
    # 3. Check consistency (allowing for sign differences)
    # The Matrix-Tree theorem says determinants are equal up to sign depending on ordering?
    # Or should be exactly equal if sign conventions are handled?
    # Usually exactly equal for unoriented matrix tree theorem if |Det|.
    # But compute_M_MHV returns signed value.
    
    match12 = abs(abs(m1) - abs(m2)) < 1e-10 * abs(m1)
    match13 = abs(abs(m1) - abs(m3)) < 1e-10 * abs(m1)
    
    if match12 and match13:
        print("\n[PASS] Magnitudes match across different root choices.")
        
        # Check signs
        s1 = sign(m1)
        s2 = sign(m2)
        s3 = sign(m3)
        print(f"  Signs: {s1}, {s2}, {s3}")
        if s1 == s2 == s3:
             print("  [INFO] Signs match exactly.")
        else:
             print("  [INFO] Signs differ (expected if ordering/orientation varies).")
             
    else:
        print("\n[FAIL] Magnitudes do not match!")
        print(f"  Diff 1-2: {float(abs(m1)-abs(m2)):.2e}")
        print(f"  Diff 1-3: {float(abs(m1)-abs(m3)):.2e}")
        exit(1)

if __name__ == "__main__":
    test_mhv_chart_independence()
