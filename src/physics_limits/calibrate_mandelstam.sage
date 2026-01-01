import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket, mandelstam_s

def calibrate_mandelstam():
    print("M1.1: Calibrating Mandelstam Normalization...")
    
    # We want to check: s_det vs sum <ij>[ji]
    # For a 2-particle channel (i, j):
    # s_ij = (p_i + p_j)^2 = 2 p_i.p_j
    # In spinor helicity: 2 p_i.p_j = <ij>[ji] = - <ij>[ij]
    
    n = 6
    lambdas, tildes = sample_spinors_from_twistor(seed=42, n=n)
    
    # 1. Check 2-particle channel
    print("\n--- 2-particle channel (0, 1) ---")
    
    # A) Matrix Determinant
    l0, lt0 = lambdas[0], tildes[0]
    l1, lt1 = lambdas[1], tildes[1]
    
    P0 = Matrix(QQ, 2, 2, [[l0[0]*lt0[0], l0[0]*lt0[1]], [l0[1]*lt0[0], l0[1]*lt0[1]]])
    P1 = Matrix(QQ, 2, 2, [[l1[0]*lt1[0], l1[0]*lt1[1]], [l1[1]*lt1[0], l1[1]*lt1[1]]])
    P_tot = P0 + P1
    
    s_det = P_tot.det()
    
    # B) Spinor Bracket via mandelstam_s
    s_mand = mandelstam_s(lambdas, tildes, 0, 1)
    
    # Correction: det(P) = <ij>[ij] = - <ij>[ji] = - s_ij
    # So s_matrix = - det(P)
    s_matrix = - s_det
    
    print(f"s_matrix (-det): {s_matrix}")
    print(f"mandelstam_s:    {s_mand}")
    
    if s_mand != 0:
        ratio = s_matrix / s_mand
        print(f"Ratio (matrix / mandelstam_s): {ratio}")
        if abs(ratio - 1) < 1e-9:
             print("  PASS: 2-particle convention matches.")
        else:
             print("  FAIL: 2-particle convention mismatch.")
        
    # 2. Check 3-particle channel (0, 1, 2)
    print("\n--- 3-particle channel (0, 1, 2) ---")
    
    l2, lt2 = lambdas[2], tildes[2]
    P2 = Matrix(QQ, 2, 2, [[l2[0]*lt2[0], l2[0]*lt2[1]], [l2[1]*lt2[0], l2[1]*lt2[1]]])
    P_tot_3 = P0 + P1 + P2
    
    s_det_3 = P_tot_3.det()
    s_matrix_3 = - s_det_3
    
    # Sum of pairs (using PHYSICAL definition <ab>[ba] = -<ab>[ab])
    s_sum_phys = QQ(0)
    indices = [0, 1, 2]
    for i in range(len(indices)):
        for j in range(i+1, len(indices)):
            a, b = indices[i], indices[j]
            term = mandelstam_s(lambdas, tildes, a, b)
            s_sum_phys += term
            
    print(f"s_matrix (-det): {s_matrix_3}")
    print(f"Sum pairs:       {s_sum_phys}")
    
    if s_sum_phys != 0:
        ratio3 = s_matrix_3 / s_sum_phys
        print(f"Ratio (matrix / SumPairs): {ratio3}")
        if abs(ratio3 - 1) < 1e-9:
             print("  PASS: 3-particle convention matches.")
        else:
             print("  FAIL: 3-particle convention mismatch.")

if __name__ == "__main__":
    calibrate_mandelstam()
