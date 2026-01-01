import sys
from sage.all import *

def bcfw_shift_spinors(lambdas, tilde_lambdas, a, b, z):
    """
    Apply BCFW shift:
    la -> la + z lb
    lb_t -> lb_t - z la_t
    """
    new_L = [list(l) for l in lambdas]
    new_Lt = [list(lt) for lt in tilde_lambdas]
    
    # la
    new_L[a][0] += z * lambdas[b][0]
    new_L[a][1] += z * lambdas[b][1]
    
    # lb_t
    new_Lt[b][0] -= z * tilde_lambdas[a][0]
    new_Lt[b][1] -= z * tilde_lambdas[a][1]
    
    # Convert back to vectors
    return [vector(QQ, l) for l in new_L], [vector(QQ, lt) for lt in new_Lt]

def solve_s_channel_bcfw(lambdas, tilde_lambdas, channel_indices, shift_a, shift_b):
    """
    Find z* such that s_channel(z*) = 0.
    """
    # Helpers needed locally or imported
    # We duplicate small helpers to avoid circular deps or complexity
    def ang_bracket(la, lb):
        return la[0]*lb[1] - la[1]*lb[0]

    def sq_bracket(lta, ltb):
        return lta[0]*ltb[1] - lta[1]*ltb[0]

    def get_s_spinor(L, Lt, i, j):
        return ang_bracket(L[i], L[j]) * sq_bracket(Lt[i], Lt[j])

    def get_s_multi(L, Lt):
        s_total = QQ(0)
        idx_list = sorted(list(channel_indices))
        for i in range(len(idx_list)):
            for j in range(i+1, len(idx_list)):
                s_total += get_s_spinor(L, Lt, idx_list[i], idx_list[j])
        return s_total

    # Evaluate at z=0 and z=1
    def get_s(z_val):
        L, Lt = bcfw_shift_spinors(lambdas, tilde_lambdas, shift_a, shift_b, z_val)
        return get_s_multi(L, Lt)

    s0 = get_s(QQ(0))
    s1 = get_s(QQ(1))
    
    slope = s1 - s0
    if slope == 0:
        return None # Constant
        
    z_star = -s0 / slope
    return z_star






