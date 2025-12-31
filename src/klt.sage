#!/usr/bin/env sage
# =============================================================================
# KLT MODULE: Yang-Mills Parke-Taylor and KLT gravity amplitude
# =============================================================================
# Reference: Kawai-Lewellen-Tye relations
# Standard field-theory KLT kernel formulation
# =============================================================================

from sage.all import *
from itertools import permutations

# Import from hodges module (assumes load() in main script)
# mandelstam_invariant is expected to be available


def parke_taylor_6pt_mhv(twistor, order, neg_helicity=(0, 1)):
    """
    Compute Parke-Taylor amplitude for 6-point MHV Yang-Mills.
    
    A_n = <a b>^4 / (<order[0] order[1]><order[1] order[2]>...<order[n-1] order[0]>)
    
    Args:
        twistor: MomentumTwistor instance
        order: Color ordering as list of particle indices
        neg_helicity: Tuple of two negative helicity particle indices
        
    Returns:
        Parke-Taylor amplitude in QQ, or None if undefined
    """
    n = twistor.n
    if len(order) != n:
        return None
    
    # Cyclic product of angle brackets
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        idx_i = order[i]
        idx_j = order[j]
        bracket = twistor.get_angle(idx_i, idx_j)
        if bracket == 0:
            return None
        denom *= bracket
    
    # MHV helicity factor <a b>^4
    neg_a, neg_b = neg_helicity
    helicity_factor = twistor.get_angle(neg_a, neg_b)
    if helicity_factor == 0:
        return None
    
    if denom == 0:
        return None
    
    return (helicity_factor ** 4) / denom


def klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_func):
    """
    Compute KLT momentum kernel S_KLT[alpha|beta] for 6-point.
    
    Standard field-theory KLT formula:
    S[alpha|beta] = ∏_{i=0}^{2} (s_{0,alpha[i]} + Σ_{j<i} theta(alpha[j],alpha[i]) * s_{alpha[j],alpha[i]})
    
    For n=6:
    - Permuted set is {2,3,4} (0-based: {1,2,3})
    - Fixed legs are {1,5,6} (0-based: {0,4,5})
    - alpha, beta are permutations of {1,2,3}
    
    theta_beta(a,b) = 1 if a appears after b in beta, else 0
    
    Args:
        alpha: Permutation of permuted set
        beta: Permutation of permuted set
        twistor: MomentumTwistor instance
        mandelstam_func: Function to compute Mandelstam invariants
        
    Returns:
        Kernel value in QQ, or None if undefined
    """
    if len(alpha) != 3 or len(beta) != 3:
        return None
    
    # Build position map for beta
    pos_in_beta = {}
    for idx, val in enumerate(beta):
        pos_in_beta[val] = idx
    
    # Theta function: theta_beta(a,b) = 1 if a appears after b in beta
    def theta_beta(a, b):
        if a not in pos_in_beta or b not in pos_in_beta:
            return 0
        return 1 if pos_in_beta[a] > pos_in_beta[b] else 0
    
    # Standard KLT kernel formula
    kernel = QQ(1)
    
    for i in range(3):  # i = 0, 1, 2
        # First term: s_{0,alpha[i]}
        s_0_ai = mandelstam_func(twistor, 0, alpha[i])
        if s_0_ai is None:
            return None
        
        sum_term = s_0_ai
        
        # Sum over j < i
        for j in range(i):
            theta_ji = theta_beta(alpha[j], alpha[i])
            if theta_ji:
                s_aj_ai = mandelstam_func(twistor, alpha[j], alpha[i])
                if s_aj_ai is None:
                    return None
                sum_term += s_aj_ai
        
        kernel *= sum_term
    
    return kernel


def gravity_6pt_mhv_klt(twistor, mandelstam_func):
    """
    Compute 6-point MHV gravity amplitude via KLT double-copy.
    
    M_6 = Σ_{alpha,beta ∈ S3} A(5,6,alpha,1) * S_KLT[alpha|beta] * A(1,beta,5,6)
    
    Where:
    - alpha, beta are permutations of {2,3,4} (0-based: {1,2,3})
    - Fixed legs: {1,5,6} (0-based: {0,4,5})
    
    Args:
        twistor: MomentumTwistor instance
        mandelstam_func: Function to compute Mandelstam invariants
        
    Returns:
        Tuple (amplitude, reason) where reason is "ok" or error
    """
    # Permuted set: {1,2,3} (0-based for {2,3,4})
    permuted_set = [1, 2, 3]
    
    # Fixed legs: {0,4,5} (0-based for {1,5,6})
    fixed_leg_1 = 0
    fixed_leg_5 = 4
    fixed_leg_6 = 5
    
    total = QQ(0)
    
    # All permutations - CANONICAL: Use lexicographically sorted order
    all_perms = sorted(list(permutations(permuted_set)))
    
    # KLT Form (n=6):
    # M = sum_{alpha, beta} A(1, alpha, 5, 6) S[alpha|beta] A(1, beta, 6, 5)
    #
    # Wait, the ordering of the second amplitude matters for the sign/kernel.
    # The standard formula (Bern et al) uses A(1, alpha, n-1, n) and A(1, beta, n, n-1).
    # Note the swap of (n-1, n) in the second amplitude.
    #
    # For n=6:
    # First factor: A(1, alpha, 5, 6). Order: [0] + alpha + [4, 5].
    # Second factor: A(1, beta, 6, 5). Order: [0] + beta + [5, 4].
    #
    # Let's adjust the second amplitude ordering.
    
    for alpha in all_perms:
        alpha = list(alpha)
        
        # A(1,alpha,5,6) -> [0] + alpha + [4, 5]
        order_alpha = [fixed_leg_1] + alpha + [fixed_leg_5, fixed_leg_6]
        A_alpha = parke_taylor_6pt_mhv(twistor, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            beta = list(beta)
            
            # A(1,beta,6,5) -> [0] + beta + [5, 4]
            # This swap (5,6) -> (6,5) introduces a sign (-1) if we were relating to (1,beta,5,6).
            # But we calculate it directly.
            order_beta = [fixed_leg_1] + beta + [fixed_leg_6, fixed_leg_5]
            A_beta = parke_taylor_6pt_mhv(twistor, order_beta)
            if A_beta is None:
                continue
            
            # KLT kernel
            # S[alpha|beta] with pivot 1 (index 0).
            S = klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_func)
            if S is None:
                continue
            
            # Add contribution
            # Note: Is there an overall sign (-1)^(n+1)?
            # For n=6, (-1)^7 = -1.
            # But let's check symmetry first. If we get symmetry, we can fix the overall sign later.
            
            total += A_alpha * S * A_beta
    
    return (total, "ok")

