import sys
import os
from sage.all import *

def ang_bracket(la, lb):
    return la[0]*lb[1] - la[1]*lb[0]

def sq_bracket(lta, ltb):
    # [i j] = u0 v1 - u1 v0
    return lta[0]*ltb[1] - lta[1]*ltb[0]

def ang_vec(a, b):
    """Angle bracket for 2-vectors: <a b> = a[0]*b[1] - a[1]*b[0]"""
    return a[0] * b[1] - a[1] * b[0]

def hodges_npt_mhv_spinor(lambdas, tilde_lambdas, neg=(0, 1), delete=(0, 1, 2)):
    """
    Generalized Hodges formula for n-point MHV gravity.
    
    Args:
        lambdas, tilde_lambdas: Lists of spinors
        neg: Indices of negative helicity particles (default (0,1))
        delete: Indices of 3 rows/cols to delete for reduced determinant (default (0,1,2))
        
    Returns:
        (Amplitude, status_string)
    """
    n = len(lambdas)
    
    # Build Phi matrix
    Phi = matrix(QQ, n, n)
    
    # Off-diagonal: Phi_{ij} = [i j] / <i j>
    for i in range(n):
        for j in range(n):
            if i != j:
                ang = ang_bracket(lambdas[i], lambdas[j])
                sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
                
                if ang == 0:
                    return None, "domain_violation_angle_bracket_offdiag"
                    
                Phi[i, j] = sq / ang
                
    # Diagonal: Reference spinor formula
    # Use generic references to avoid accidental zeros
    lx = vector(QQ, [1, 2])
    ly = vector(QQ, [3, 1])
    
    # Check if valid references (not orthogonal to any lambda)
    for i in range(n):
        if ang_bracket(lambdas[i], lx) == 0:
            lx = vector(QQ, [1, 1]) 
        if ang_bracket(lambdas[i], ly) == 0:
            ly = vector(QQ, [1, -1])

    for i in range(n):
        ang_ix = ang_bracket(lambdas[i], lx)
        ang_iy = ang_bracket(lambdas[i], ly)
        
        if ang_ix == 0 or ang_iy == 0:
             return None, "domain_violation_reference_spinor"
             
        diag_sum = QQ(0)
        for j in range(n):
            if j == i: continue
            
            ang_jx = ang_bracket(lambdas[j], lx)
            ang_jy = ang_bracket(lambdas[j], ly)
            
            term = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= term
            
        Phi[i, i] = diag_sum
        
    # Reduced determinant
    rows_to_delete = list(delete)
    cols_to_delete = rows_to_delete
    
    rows_keep = [r for r in range(n) if r not in rows_to_delete]
    cols_keep = [c for c in range(n) if c not in cols_to_delete]
    
    Phi_red = Phi[rows_keep, cols_keep]
    det_Phi_red = Phi_red.det()
    
    # Normalization factor for deleted set {r1, r2, r3}
    # norm = (<r1 r2><r2 r3><r3 r1>)^2
    r1, r2, r3 = rows_to_delete
    ang_12 = ang_bracket(lambdas[r1], lambdas[r2])
    ang_23 = ang_bracket(lambdas[r2], lambdas[r3])
    ang_31 = ang_bracket(lambdas[r3], lambdas[r1])
    
    norm_factor = (ang_12 * ang_23 * ang_31)**2
    if norm_factor == 0:
        return None, "domain_violation_norm_factor"
        
    det_prime = det_Phi_red / norm_factor
    
    # Helicity factor <a b>^8
    a, b = neg
    h_factor = ang_bracket(lambdas[a], lambdas[b])**8
    
    # Note: There might be a global sign depending on n and deletion set.
    # Hodges paper suggests sign is related to permutation of deletion set.
    # For symmetric deletion (0,1,2), the sign convention usually matches.
    # We return the raw value here.
    
    return det_prime * h_factor, "ok"

def hodges_6pt_mhv_spinor(lambdas, tilde_lambdas, deletion_set=None):
    """Legacy wrapper for 6pt."""
    if deletion_set is None:
        deletion_set = (0, 1, 2)
    return hodges_npt_mhv_spinor(lambdas, tilde_lambdas, neg=(0, 1), delete=deletion_set)

