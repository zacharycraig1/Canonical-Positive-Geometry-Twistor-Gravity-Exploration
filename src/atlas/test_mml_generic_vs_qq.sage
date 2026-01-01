import sys
import os
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
from src.posgeom.moment_map_laplacian import MomentMapLaplacian
from src.posgeom.forest_polytope import get_forest_exponents

def test_mml_generic_vs_qq():
    print("Running Unit Test: MML Generic vs QQ...")
    
    n = 6
    roots = [0, 1, 2]
    
    # 1. Setup MML and Evaluator
    evaluator = FiberGaugeEvaluator(n)
    # We can access _compute_X_H_generic directly from evaluator instance
    
    # Create MML instance directly to compare against
    exponents, edge_order = get_forest_exponents(n, roots)
    mml = MomentMapLaplacian(n, roots, edge_order)
    
    # 2. Generate integer inputs
    # MML expects z_values corresponding to edges
    num_edges = len(edge_order)
    z_vals = [QQ(rnd.randint(1, 10)) for _ in range(num_edges)]
    
    print(f"Testing with z_vals: {z_vals}")
    
    # 3. Compute using original MML (QQ)
    print("Computing with MML.compute_X_H (QQ)...")
    X1, H1 = mml.compute_X_H(z_vals)
    
    # 4. Compute using Generic Evaluator (QQ)
    print("Computing with Evaluator._compute_X_H_generic (QQ)...")
    X2, H2 = evaluator._compute_X_H_generic(mml, z_vals, QQ)
    
    # 5. Compare
    print("\nComparing Results...")
    
    # X check
    diff_X = vector(QQ, X1) - vector(QQ, X2)
    if diff_X.norm() == 0:
        print("  [PASS] X vectors match exactly.")
    else:
        print("  [FAIL] X vectors differ!")
        print(f"  X1: {X1}")
        print(f"  X2: {X2}")
        print(f"  Diff: {diff_X}")
        exit(1)
        
    # H check
    diff_H = matrix(QQ, H1) - matrix(QQ, H2)
    if diff_H.norm() == 0:
        print("  [PASS] H matrices match exactly.")
    else:
        print("  [FAIL] H matrices differ!")
        print(f"  H1 norm: {matrix(QQ, H1).norm()}")
        print(f"  H2 norm: {matrix(QQ, H2).norm()}")
        print(f"  Diff norm: {diff_H.norm()}")
        exit(1)

    print("\nSuccess! Generic implementation matches reference on QQ.")

if __name__ == "__main__":
    test_mml_generic_vs_qq()
