import sys
import os
from sage.all import *

def ang_bracket(la, lb):
    return la[0]*lb[1] - la[1]*lb[0]

def sq_bracket(lta, ltb):
    return lta[0]*ltb[1] - lta[1]*ltb[0]

def hodges_reduced_det(Phi, lambdas, rows_to_delete):
    """
    Compute reduced determinant of Phi with canonical normalization.
    
    normalization = (<r1 r2><r2 r3><r3 r1>)^2
    where {r1, r2, r3} are the deleted rows.
    
    The sign of the reduced determinant depends on the permutation of the deleted rows relative to (0, 1, 2).
    However, if we fix the order of rows_to_delete to be sorted, the sign is well-defined.
    Wait, Hodges formula: det_red = (-1)^{i+j+...} M_{ijk}^{ijk} / (<ij><jk><ki>)^2?
    The standard formula:
    M_n = (-1)^{n+1} \sigma_{ijk} \det(\Phi^{ijk}_{ijk})
    where \sigma_{ijk} = 1 / (<ij><jk><ki>)^2.
    
    Wait, (-1)^{n+1} is global.
    We just want independence from deletion set.
    """
    n = Phi.nrows()
    rows_to_delete = sorted(list(rows_to_delete))
    cols_to_delete = rows_to_delete
    
    rows_keep = [r for r in range(n) if r not in rows_to_delete]
    cols_keep = [c for c in range(n) if c not in cols_to_delete]
    
    Phi_red = Phi[rows_keep, cols_keep]
    det_Phi_red = Phi_red.det()
    
    r1, r2, r3 = rows_to_delete
    ang_12 = ang_bracket(lambdas[r1], lambdas[r2])
    ang_23 = ang_bracket(lambdas[r2], lambdas[r3])
    ang_31 = ang_bracket(lambdas[r3], lambdas[r1])
    
    norm_factor = (ang_12 * ang_23 * ang_31)**2
    
    if norm_factor == 0:
        return None, "domain_violation_norm_factor"
        
    det_prime = det_Phi_red / norm_factor
    
    return det_prime, "ok"

def hodges_npt_mhv_canonical(lambdas, tilde_lambdas, negative_indices, deletion_set=None):
    """
    Compute MHV amplitude using Hodges formula, checking deletion set independence.
    """
    n = len(lambdas)
    
    # Build Phi
    Phi = matrix(QQ, n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                ang = ang_bracket(lambdas[i], lambdas[j])
                sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
                if ang == 0: return None, "domain_violation_angle_bracket_offdiag"
                Phi[i, j] = sq / ang
                
    # Diagonal
    lx = vector(QQ, [1, 2])
    ly = vector(QQ, [3, 1])
    
    for i in range(n):
        if ang_bracket(lambdas[i], lx) == 0: lx = vector(QQ, [1, 1])
        if ang_bracket(lambdas[i], ly) == 0: ly = vector(QQ, [1, -1])
    
    for i in range(n):
        diag_sum = QQ(0)
        ang_ix = ang_bracket(lambdas[i], lx)
        ang_iy = ang_bracket(lambdas[i], ly)
        if ang_ix == 0 or ang_iy == 0: return None, "domain_violation_ref"
        
        for j in range(n):
            if j == i: continue
            ang_jx = ang_bracket(lambdas[j], lx)
            ang_jy = ang_bracket(lambdas[j], ly)
            term = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= term
        Phi[i, i] = diag_sum
        
    # Compute reduced det
    if deletion_set is None:
        target_set = (0, 1, 2)
    else:
        target_set = deletion_set
        
    val_main, status = hodges_reduced_det(Phi, lambdas, target_set)
    
    if status != "ok":
        # Only fallback if user didn't specify a set (if they did, they want THAT one)
        if deletion_set is None:
             val_main, status = hodges_reduced_det(Phi, lambdas, (n-3, n-2, n-1))
             if status != "ok":
                 return None, status
        else:
             return None, status

    # Helicity factor <a b>^8
    a, b = negative_indices
    h_factor = ang_bracket(lambdas[a], lambdas[b])**8
    
    return val_main * h_factor, "ok"

