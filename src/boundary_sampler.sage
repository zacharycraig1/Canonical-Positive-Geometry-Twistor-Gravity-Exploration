from sage.all import *

try:
    from src.spinor_sampling import sample_spinor_helicity_conserving
    from src.kinematics_map import SpinorHelicityAdapter, ang_bracket, sq_bracket
except ImportError:
    # Assuming loaded via load()
    pass

def solve_boundary_shift(lambdas, tilde_lambdas, S, shift_pair):
    """
    Find shift parameter z such that s_S(z) = 0.
    Shift is <a, b]:
      tilde_lambda_a -> tilde_lambda_a - z * tilde_lambda_b
      lambda_b -> lambda_b + z * lambda_a
      
    Args:
        lambdas, tilde_lambdas: Initial spinors
        S: Set of indices defining the channel
        shift_pair: Tuple (a, b) where a in S, b not in S (or vice versa)
        
    Returns:
        z_pole (QQ)
    """
    a, b = shift_pair
    
    # Ensure correct direction: a in S, b not in S
    # If a not in S and b in S, the momentum P_S(z) shifts by +k instead of -k
    # P_S(z) = sum_{i in S} p_i(z)
    # If a in S: p_a(z) = p_a - z |a><b|. Contribution -z|a><b|.
    # If b in S: p_b(z) = p_b + z |a><b|. Contribution +z|a><b|.
    
    sign = 0
    if a in S and b not in S:
        # P_S(z) = P_S - z |a><b|
        # s_S(z) = s_S + z * matrix_elem (due to sq convention)
        # z = -s_S / matrix_elem
        # denom needs to be matrix_elem, so sign=1
        sign = 1
    elif b in S and a not in S:
        # P_S(z) = P_S + z |a><b|
        # s_S(z) = s_S - z * matrix_elem
        # z = s_S / matrix_elem
        # denom needs to be -matrix_elem, so sign=-1
        sign = -1
    else:
        raise ValueError("Shift pair must bridge S and complement")

    # Current s_S = P_S^2
    # P_S = sum_{i in S} p_i
    # We don't need full vector, just s_S and <a|P_S|b]
    
    # Compute s_S(0)
    # s_S = sum_{i<j in S} s_{ij}
    s_S_0 = QQ(0)
    S_list = sorted(list(S))
    for idx1 in range(len(S_list)):
        for idx2 in range(idx1 + 1, len(S_list)):
            i, j = S_list[idx1], S_list[idx2]
            # s_ij = <ij>[ij]
            # We assume [ij] = det(tilde_i, tilde_j)
            s_ij = ang_bracket(lambdas, i, j) * sq_bracket(tilde_lambdas, i, j)
            s_S_0 += s_ij
            
    # Compute <a|P_S|b] = sum_{i in S} <a i> [i b]
    matrix_elem = QQ(0)
    for i in S:
        val = ang_bracket(lambdas, a, i) * sq_bracket(tilde_lambdas, i, b)
        matrix_elem += val
        
    # Equation: s_S(z) = (P + sign*z*k)^2 = P^2 + 2*sign*z*(P.k)
    # k = |a>[b| = lambda_a tilde_lambda_b
    # 2 P.k = 2 P_{mu} 0.5 <a|sigma_mu|b] ? No.
    # 2 P.k = <a|P|b]
    # So s_S(z) = s_S(0) + sign * z * <a|P_S|b]
    
    denom = sign * matrix_elem
    
    if denom == 0:
        return None # Cannot shift to boundary with this pair
        
    z_pole = -s_S_0 / denom
    return z_pole

def apply_shift(lambdas, tilde_lambdas, a, b, z):
    """
    Apply BCFW shift <a, b] by amount z.
    Returns new lists of spinors.
    """
    new_lambdas = [v for v in lambdas]
    new_tilde = [v for v in tilde_lambdas]
    
    # tilde_lambda_a -> tilde_lambda_a - z * tilde_lambda_b
    new_tilde[a] = tilde_lambdas[a] - z * tilde_lambdas[b]
    
    # lambda_b -> lambda_b + z * lambda_a
    new_lambdas[b] = lambdas[b] + z * lambdas[a]
    
    return new_lambdas, new_tilde

def sample_near_boundary(S, epsilon, seed=None):
    """
    Sample kinematics near boundary s_S = 0.
    
    Args:
        S: Tuple/list of indices
        epsilon: Distance from boundary (s_S = epsilon, approx?)
                 Actually we set z = z_pole + epsilon.
                 s_S(z_pole+eps) ~ eps * (slope).
                 So s_S will be proportional to epsilon.
        seed: Random seed
        
    Returns:
        (lambdas, tilde_lambdas)
    """
    import random
    if seed is not None:
        random.seed(seed)
        
    # 1. Generic sample
    # Use a seed for the underlying sampler derived from main seed
    gen_seed = random.randint(1, 1000000)
    res = sample_spinor_helicity_conserving(n=6, seed=gen_seed)
    if res is None: return None
    lambdas, tilde_lambdas = res
    
    # 2. Pick shift pair
    # a in S, b not in S
    n = 6
    S_set = set(S)
    complement = [x for x in range(n) if x not in S_set]
    
    if not complement:
        return None # S is full set? s_S is always 0 by conservation?
        
    a = random.choice(list(S_set))
    b = random.choice(complement)
    
    # 3. Solve for pole
    z_pole = solve_boundary_shift(lambdas, tilde_lambdas, S_set, (a, b))
    
    if z_pole is None:
        return None
        
    # 4. Apply shift z = z_pole + epsilon
    # We handle epsilon as a generic parameter?
    # No, we want concrete numbers.
    # If epsilon is float, result is float. If Rational, result is Rational.
    z_target = z_pole + epsilon
    
    return apply_shift(lambdas, tilde_lambdas, a, b, z_target)


