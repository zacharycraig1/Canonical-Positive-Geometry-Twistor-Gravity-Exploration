import sys
from sage.all import *
from itertools import permutations

def parke_taylor_n_pt(twistor, order, neg_helicity=(0, 1)):
    n = len(order)
    if n != twistor.n:
        return None
        
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        idx_i = order[i]
        idx_j = order[j]
        bracket = twistor.get_angle(idx_i, idx_j)
        if bracket == 0: return None
        denom *= bracket
        
    neg_a, neg_b = neg_helicity
    h_factor = twistor.get_angle(neg_a, neg_b)
    if h_factor == 0: return None
    
    return (h_factor ** 4) / denom

def parke_taylor_6pt_mhv(twistor, order, neg_helicity=(0, 1)):
    """
    Compute Parke-Taylor amplitude for 6-point MHV Yang-Mills.
    
    A_n = <a b>^4 / (<order[0] order[1]><order[1] order[2]>...<order[n-1] order[0]>)
    """
    return parke_taylor_n_pt(twistor, order, neg_helicity)


def klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_func):
    """
    Compute KLT momentum kernel S_KLT[alpha|beta] for 6-point.
    
    Standard field-theory KLT formula:
    S[alpha|beta] = ∏_{i=0}^{2} (s_{0,alpha[i]} + Σ_{j<i} theta(alpha[j],alpha[i]) * s_{alpha[j],alpha[i]})
    """
    if len(alpha) != 3 or len(beta) != 3:
        return None
    
    # Build position map for beta
    pos_in_beta = {}
    for idx, val in enumerate(beta):
        pos_in_beta[val] = idx
    
    # Theta function: theta_beta(a,b) = 1 if a appears BEFORE b in beta
    # Changed from > to < to match n=5 check
    def theta_beta(a, b):
        if a not in pos_in_beta or b not in pos_in_beta:
            return 0
        return 1 if pos_in_beta[a] < pos_in_beta[b] else 0
    
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
            S = klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_func)
            if S is None:
                continue
            
            total += A_alpha * S * A_beta
    
    return (total, "ok")

