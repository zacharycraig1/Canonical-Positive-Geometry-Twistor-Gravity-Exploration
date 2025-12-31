def sample_spinor_helicity_conserving(n=6, seed=None):
    """
    Sample spinor helicity variables (lambda, tilde_lambda) satisfying momentum conservation.
    
    Args:
        n: Number of particles (default 6)
        seed: Random seed
        
    Returns:
        (lambdas, tilde_lambdas)
        lambdas: list of n 2-vectors (QQ)
        tilde_lambdas: list of n 2-vectors (QQ)
        
    Momentum conservation: sum_i lambda_i * tilde_lambda_i^T = 0
    """
    import random
    if seed is not None:
        random.seed(int(seed))
        
    # 1. Sample random lambdas (generic)
    lambdas = []
    for i in range(n):
        l = vector(QQ, [QQ(random.randint(-100, 100)), QQ(random.randint(-100, 100))])
        # Ensure no lambda is zero and generic enough (no collinear pairs ideally)
        while l == 0:
            l = vector(QQ, [QQ(random.randint(-100, 100)), QQ(random.randint(-100, 100))])
        lambdas.append(l)
        
    # 2. Sample random tilde_lambdas for i=0..n-3
    tilde_lambdas = [None] * n
    for i in range(n-2):
        l_tilde = vector(QQ, [QQ(random.randint(-100, 100)), QQ(random.randint(-100, 100))])
        tilde_lambdas[i] = l_tilde
        
    # 3. Solve for tilde_lambda_{n-2} and tilde_lambda_{n-1} (indices n-2, n-1)
    # Equation: sum_{i=0}^{n-1} lambda_i * tilde_lambda_i^T = 0
    # Let L = [lambda_0 ... lambda_{n-1}] (2xn), Lt = [tilde_lambda_0 ...]^T (nx2)
    # L * Lt = 0 (2x2 zero matrix)
    
    # Rearrange: sum_{i=0}^{n-3} lambda_i * tilde_lambda_i^T + lambda_{n-2}*tilde_{n-2}^T + lambda_{n-1}*tilde_{n-1}^T = 0
    # RHS = - sum_{i=0}^{n-3} ...
    
    P_known = matrix(QQ, 2, 2, 0)
    for i in range(n-2):
        P_known += matrix(QQ, 2, 2, [lambdas[i][0]*tilde_lambdas[i][0], lambdas[i][0]*tilde_lambdas[i][1],
                                     lambdas[i][1]*tilde_lambdas[i][0], lambdas[i][1]*tilde_lambdas[i][1]])
                                     
    RHS = -P_known
    
    # We need to solve:
    # lambda_{n-2} * tilde_{n-2}^T + lambda_{n-1} * tilde_{n-1}^T = RHS
    # This is a linear system for the 4 components of tilde_{n-2}, tilde_{n-1}
    # It decomposes into 2 separate systems for the two columns of tilde_lambdas.
    
    # Let tilde_lambda_k = (u_k, v_k).
    # Col 1 of RHS (RHS_00, RHS_10) comes from u's.
    # Col 2 of RHS (RHS_01, RHS_11) comes from v's.
    
    # System for u_{n-2}, u_{n-1}:
    # lambda_{n-2}[0]*u_{n-2} + lambda_{n-1}[0]*u_{n-1} = RHS[0,0]
    # lambda_{n-2}[1]*u_{n-2} + lambda_{n-1}[1]*u_{n-1} = RHS[1,0]
    
    # Matrix M = [lambda_{n-2}, lambda_{n-1}] (columns)
    M = matrix(QQ, 2, 2, [lambdas[n-2][0], lambdas[n-1][0], 
                          lambdas[n-2][1], lambdas[n-1][1]])
                          
    if M.det() == 0:
        return None # Failed, last two lambdas collinear
        
    # Solve for u's
    b_u = vector(QQ, [RHS[0,0], RHS[1,0]])
    u_sol = M.solve_right(b_u)
    
    # Solve for v's
    b_v = vector(QQ, [RHS[0,1], RHS[1,1]])
    v_sol = M.solve_right(b_v)
    
    tilde_lambdas[n-2] = vector(QQ, [u_sol[0], v_sol[0]])
    tilde_lambdas[n-1] = vector(QQ, [u_sol[1], v_sol[1]])
    
    return lambdas, tilde_lambdas

