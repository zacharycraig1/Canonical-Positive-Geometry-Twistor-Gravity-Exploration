#!/usr/bin/env sage
"""
Pushforward n=6 Failure Analysis
================================

This script systematically analyzes why the saddle pushforward
works for n=4,5 but fails for n=6.

Hypotheses:
1. Wrong normalization (constant scaling factor)
2. Missing Jacobian determinant
3. Incorrect saddle counting
4. Wrong moment map formulation
"""

from sage.all import *
import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.posgeom.forest_polytope import get_forest_polynomial, get_forest_exponents
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket

# Try to import the saddle pushforward (numpy-based)
try:
    from src.posgeom.saddle_pushforward import compute_pushforward_saddle, moment_map_and_jacobian
    HAS_PUSHFORWARD = True
except ImportError:
    HAS_PUSHFORWARD = False
    print("Warning: saddle_pushforward not available")


def get_reference_spinors():
    """Standard reference spinors."""
    return vector(QQ, [1, 0]), vector(QQ, [0, 1])


def compute_forest_polynomial_value(n, roots, lambdas, tildes, x, y):
    """
    Evaluate F_{n,R}(z) at the given kinematics.
    """
    # Get the forest polynomial
    F_poly = get_forest_polynomial(n, roots)
    
    # Compute z_ij values
    z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
    
    # Substitute into polynomial
    R_ring = F_poly.parent()
    eval_dict = {}
    for var_name in R_ring.variable_names():
        if var_name in z_map:
            eval_dict[R_ring(var_name)] = z_map[var_name]
    
    F_val = F_poly.subs(eval_dict)
    
    if hasattr(F_val, 'is_constant') and F_val.is_constant():
        return F_val.constant_coefficient()
    return F_val


def test_ratio_constancy(n_trials=50):
    """
    Test if pushforward / Hodges ratio is constant across kinematic points.
    
    If CONSTANT → we just need to find the normalization factor
    If VARYING → the map itself is wrong
    """
    print("="*70)
    print("HYPOTHESIS 1: Is the ratio constant?")
    print("="*70)
    
    n = 6
    roots = [0, 1, 2]
    x, y = get_reference_spinors()
    
    # Cache the polynomial and exponents
    print("Building forest polynomial for n=6...")
    F_poly = get_forest_polynomial(n, roots)
    exponents, edge_order = get_forest_exponents(n, roots)
    
    # Convert to numpy for pushforward
    poly_exponents = np.array(exponents)
    poly_coeffs = np.ones(len(exponents))  # All coefficients are 1
    
    print(f"Forest polynomial has {len(exponents)} terms")
    print(f"Dimension of parameter space: {len(edge_order)}")
    print()
    
    ratios = []
    hodges_vals = []
    forest_vals = []
    
    for trial in range(n_trials):
        try:
            # Sample kinematics
            lambdas, tildes = sample_spinors_from_twistor(n=n, seed=100+trial)
            
            # Compute z_ij values
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
            
            # Get z values in canonical edge order
            z_vec = []
            for edge in edge_order:
                i, j = edge
                key = f"z_{i}_{j}"
                z_vec.append(float(z_map[key]))
            z_vec = np.array(z_vec)
            
            # Skip if any z is non-positive (outside domain)
            if np.any(z_vec <= 0):
                continue
            
            log_z = np.log(z_vec)
            
            # Compute moment map X = grad_logz log(F)
            X, J = moment_map_and_jacobian(log_z, poly_coeffs, poly_exponents)
            det_J = np.linalg.det(J)
            
            if abs(det_J) < 1e-15:
                continue
            
            pushforward_val = 1.0 / det_J
            
            # Compute Hodges reference
            M_hodges, status = reconstruct_mhv_from_laplacian(
                lambdas, tildes, x, y, roots=tuple(roots)
            )
            
            if M_hodges is None or M_hodges == 0:
                continue
            
            # Also compute forest polynomial value for comparison
            F_val = compute_forest_polynomial_value(n, roots, lambdas, tildes, x, y)
            
            ratio = pushforward_val / float(M_hodges)
            ratios.append(ratio)
            hodges_vals.append(float(M_hodges))
            forest_vals.append(float(F_val))
            
        except Exception as e:
            continue
    
    print(f"Collected {len(ratios)} valid data points")
    print()
    
    if len(ratios) < 5:
        print("ERROR: Not enough valid data points")
        return None
    
    ratios = np.array(ratios)
    
    # Statistics
    mean_ratio = np.mean(ratios)
    std_ratio = np.std(ratios)
    min_ratio = np.min(ratios)
    max_ratio = np.max(ratios)
    
    print(f"Ratio statistics:")
    print(f"  Mean:   {mean_ratio:.6e}")
    print(f"  Std:    {std_ratio:.6e}")
    print(f"  Min:    {min_ratio:.6e}")
    print(f"  Max:    {max_ratio:.6e}")
    print(f"  CV:     {std_ratio/abs(mean_ratio):.6e}")  # Coefficient of variation
    print()
    
    # Decision
    cv = std_ratio / abs(mean_ratio) if mean_ratio != 0 else float('inf')
    
    if cv < 0.01:  # Less than 1% variation
        print("[RESULT] Ratio is CONSTANT!")
        print(f"  Missing normalization factor: {1/mean_ratio:.6e}")
        print("  → Pushforward IS correct up to normalization")
        return mean_ratio
    else:
        print("[RESULT] Ratio is NOT constant")
        print("  → The moment map or pushforward is fundamentally wrong")
        return None


def analyze_jacobian_properties(n_trials=20):
    """
    Analyze properties of the Jacobian matrix at various kinematic points.
    """
    print()
    print("="*70)
    print("HYPOTHESIS 2: Jacobian Properties")
    print("="*70)
    
    n = 6
    roots = [0, 1, 2]
    x, y = get_reference_spinors()
    
    exponents, edge_order = get_forest_exponents(n, roots)
    poly_exponents = np.array(exponents)
    poly_coeffs = np.ones(len(exponents))
    
    print(f"Jacobian dimension: {len(edge_order)} x {len(edge_order)}")
    print()
    
    for trial in range(n_trials):
        try:
            lambdas, tildes = sample_spinors_from_twistor(n=n, seed=200+trial)
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
            
            z_vec = []
            for edge in edge_order:
                i, j = edge
                z_vec.append(float(z_map[f"z_{i}_{j}"]))
            z_vec = np.array(z_vec)
            
            if np.any(z_vec <= 0):
                continue
            
            log_z = np.log(z_vec)
            X, J = moment_map_and_jacobian(log_z, poly_coeffs, poly_exponents)
            
            # Eigenvalue analysis
            eigenvalues = np.linalg.eigvals(J)
            eigenvalues_real = np.real(eigenvalues)
            
            det_J = np.linalg.det(J)
            cond = np.linalg.cond(J)
            
            print(f"Trial {trial}:")
            print(f"  det(J) = {det_J:.6e}")
            print(f"  cond(J) = {cond:.6e}")
            print(f"  min(λ) = {np.min(eigenvalues_real):.6e}")
            print(f"  max(λ) = {np.max(eigenvalues_real):.6e}")
            
            # Check if J is positive definite (should be for covariance matrix)
            if np.all(eigenvalues_real > 0):
                print(f"  J is POSITIVE DEFINITE ✓")
            else:
                print(f"  J is NOT positive definite ✗")
            print()
            
        except Exception as e:
            print(f"Trial {trial}: Error - {e}")
            continue


def compare_n4_n5_n6():
    """
    Compare pushforward behavior across n=4, 5, 6.
    """
    print()
    print("="*70)
    print("COMPARISON: n=4, n=5, n=6")
    print("="*70)
    
    for n in [4, 5, 6]:
        roots = list(range(min(3, n)))  # [0,1,2] or [0,1] for n=4
        
        print(f"\nn={n}, roots={roots}")
        print("-"*40)
        
        try:
            exponents, edge_order = get_forest_exponents(n, roots)
            print(f"  Polynomial terms: {len(exponents)}")
            print(f"  Edge dimensions: {len(edge_order)}")
            
            # Expected: Newton polytope dimension
            # For forest polytope of K_n, dimension is roughly (n choose 2) - n + 1
            expected_dim = len(edge_order)
            print(f"  Newton polytope intrinsic dim: ≤ {expected_dim}")
            
        except Exception as e:
            print(f"  Error: {e}")


def main():
    print("="*70)
    print("PUSHFORWARD N=6 FAILURE ANALYSIS")
    print("="*70)
    print()
    
    if not HAS_PUSHFORWARD:
        print("ERROR: saddle_pushforward module not available")
        print("Cannot run numerical analysis")
        return
    
    # Test 1: Is ratio constant?
    ratio = test_ratio_constancy(n_trials=30)
    
    # Test 2: Jacobian properties
    analyze_jacobian_properties(n_trials=10)
    
    # Test 3: Compare n=4,5,6
    compare_n4_n5_n6()
    
    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    
    if ratio is not None:
        print(f"✓ Ratio is constant: {ratio:.6e}")
        print(f"  → Normalization factor needed: {1/ratio:.6e}")
        print("  → Proceed to investigate normalization origin")
    else:
        print("✗ Ratio is not constant or data insufficient")
        print("  → Pushforward approach may be fundamentally limited")
        print("  → Recommend proceeding to BCFW Amplituhedron (Option A)")


if __name__ == '__main__':
    main()


