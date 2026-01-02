#!/usr/bin/env sage
# =============================================================================
# HODGES SIGMA SIGN COMPUTATION
# =============================================================================
# Computes σ(ijk,rst) - the sign of permutation mapping row ordering to column ordering
# =============================================================================

from sage.all import *
from sage.combinat.permutation import Permutation

def hodges_sigma(I, J, n):
    """
    Compute Hodges sigma sign σ(ijk,rst).
    
    Args:
        I: List of deleted row indices [i,j,k] (0-indexed)
        J: List of deleted column indices [r,s,t] (0-indexed)
        n: Total number of particles (6 for n=6)
        
    Returns:
        Sign: +1 or -1
    """
    # Convert to 1-indexed for the algorithm (Hodges uses 1-indexed)
    I_1based = [i + 1 for i in I]
    J_1based = [j + 1 for j in J]
    
    # Get remaining labels (1-indexed)
    all_labels = list(range(1, n + 1))
    restI = sorted([x for x in all_labels if x not in I_1based])
    restJ = sorted([x for x in all_labels if x not in J_1based])
    
    # Build ordering lists
    A = I_1based + restI
    B = J_1based + restJ
    
    # Compute permutation that maps A to B
    # p[m] = position of A[m] in B
    posB = {val: idx for idx, val in enumerate(B)}
    p_list = [posB[val] for val in A]
    
    # Convert to 0-indexed for Permutation (Sage uses 0-indexed)
    p_0based = [x for x in p_list]
    
    # Compute sign using Sage's Permutation
    try:
        perm = Permutation(p_0based)
        return perm.sign()
    except:
        # Fallback: compute sign by inversion count
        inversions = 0
        for i in range(len(p_0based)):
            for j in range(i + 1, len(p_0based)):
                if p_0based[i] > p_0based[j]:
                    inversions += 1
        return 1 if inversions % 2 == 0 else -1










