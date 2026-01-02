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

def mandelstam_s(lambdas, tilde_lambdas, i, j):
    """
    Compute Mandelstam invariant s_ij with fixed convention.
    s_ij = (p_i + p_j)^2 = 2 p_i.p_j
    Convention: s_ij = <ij>[ji] = - <ij>[ij]
    """
    # <ij>[ji]
    return ang_bracket(lambdas[i], lambdas[j]) * sq_bracket(tilde_lambdas[j], tilde_lambdas[i])

try:
    from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
except ImportError:
    from chy_oracle.hodges_reduced import hodges_npt_mhv_canonical

def hodges_npt_mhv_spinor(lambdas, tilde_lambdas, neg=(0, 1), delete=(0, 1, 2)):
    """
    Generalized Hodges formula for n-point MHV gravity.
    Now delegates to hodges_npt_mhv_canonical to ensure consistency.
    The 'delete' argument is ignored (handled internally by canonical).
    """
    return hodges_npt_mhv_canonical(lambdas, tilde_lambdas, neg)

def hodges_6pt_mhv_spinor(lambdas, tilde_lambdas, deletion_set=None):
    """Legacy wrapper for 6pt."""
    if deletion_set is None:
        deletion_set = (0, 1, 2)
    return hodges_npt_mhv_spinor(lambdas, tilde_lambdas, neg=(0, 1), delete=deletion_set)

