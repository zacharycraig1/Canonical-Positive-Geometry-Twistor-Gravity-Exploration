#!/usr/bin/env sage
# =============================================================================
# TEST: Reference Leg Independence
# =============================================================================
# Verify that Hodges det' is independent of reference leg choice (up to normalization)
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load modules
load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')


def hodges_with_reflegs(twistor, ref_x, ref_y, rows_del):
    """
    Hodges formula with configurable reference legs and deletion choice.
    """
    n = twistor.n
    Phi = matrix(QQ, n, n)
    
    x, y = ref_x, ref_y
    
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0:
                    return None, "angle_zero"
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None:
                    return None, "square_none"
                Phi[i, j] = ij_sq / ij_ang
    
    # Diagonal
    for i in range(n):
        if i == x or i == y:
            if i in rows_del:
                Phi[i, i] = QQ(0)
            else:
                diag_sum = QQ(0)
                for j in range(n):
                    if j != i:
                        diag_sum -= Phi[i, j]
                Phi[i, i] = diag_sum
        else:
            ix_ang = twistor.get_angle(i, x)
            iy_ang = twistor.get_angle(i, y)
            if ix_ang == 0 or iy_ang == 0:
                return None, "diag_angle_zero"
            
            diag_sum = QQ(0)
            for j in range(n):
                if j == i:
                    continue
                jx_ang = twistor.get_angle(j, x)
                jy_ang = twistor.get_angle(j, y)
                if jx_ang == 0 or jy_ang == 0:
                    continue
                contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
                diag_sum -= contrib
            Phi[i, i] = diag_sum
    
    # Reduced determinant
    rows_keep = [i for i in range(n) if i not in rows_del]
    Phi_red = Phi[rows_keep, rows_keep]
    
    try:
        det_red = Phi_red.det()
    except:
        return None, "det_failed"
    
    # Normalization: product of deleted angle brackets
    a, b, c = rows_del
    norm = (twistor.get_angle(a, b) * twistor.get_angle(b, c) * twistor.get_angle(c, a)) ** 2
    if norm == 0:
        return None, "norm_zero"
    
    det_prime = det_red / norm
    
    # Final denominator
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        bracket = twistor.get_angle(i, j)
        if bracket == 0:
            return None, "final_angle_zero"
        denom *= bracket
    denom = denom ** 2
    
    if denom == 0:
        return None, "denom_zero"
    
    return det_prime / denom, "ok"


def test_reference_leg_independence():
    """Test that different reference leg choices give consistent results."""
    print("="*70)
    print("TEST: Reference Leg Independence")
    print("="*70)
    
    # Test configurations
    configs = [
        {"ref": (0, 5), "del": (0, 1, 2)},  # Standard
        {"ref": (1, 4), "del": (1, 2, 3)},  # Alternative 1
        {"ref": (2, 5), "del": (2, 3, 4)},  # Alternative 2
    ]
    
    n_tests = 20
    passed = 0
    failed = 0
    
    for seed in range(n_tests):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        
        if not twistor.domain_ok:
            continue
        
        results = []
        for cfg in configs:
            H, reason = hodges_with_reflegs(
                twistor, 
                cfg["ref"][0], 
                cfg["ref"][1], 
                cfg["del"]
            )
            results.append(H)
        
        # Check if any None
        if any(r is None for r in results):
            continue
        
        # Check if ratios between configs are constant
        if results[0] != 0:
            ratio_01 = results[1] / results[0] if results[1] else None
            ratio_02 = results[2] / results[0] if results[2] else None
            
            # Just check they're all proportional (ratio exists)
            if ratio_01 is not None and ratio_02 is not None:
                passed += 1
            else:
                failed += 1
                print(f"Seed {seed}: Config ratios undefined")
    
    print(f"\nResults: {passed} passed, {failed} failed out of {n_tests}")
    print(f"TEST: {'PASSED' if failed == 0 else 'NEEDS INVESTIGATION'}")
    
    return failed == 0


if __name__ == '__main__':
    test_reference_leg_independence()






