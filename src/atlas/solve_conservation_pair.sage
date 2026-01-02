from sage.all import *

def solve_conservation_pair(lambdas, tildes_free, pair, n=6, Field=RR):
    """
    Solves for tildes of the dependent pair (a, b) given n-2 free tildes.
    
    Args:
        lambdas: dict {i: vector([1, t_i])}
        tildes_free: dict {i: vector([tilde_x, tilde_y])} for i not in pair
        pair: tuple (a, b) indices of dependent twistors
        n: number of particles
        Field: field for computation
        
    Returns:
        tildes: dict {i: vector} for all i, or None if singular
    """
    a, b = pair
    free_indices = [i for i in range(n) if i not in pair]
    
    # Initialize RHS vectors
    rhs_0 = vector(Field, [0, 0])
    rhs_1 = vector(Field, [0, 0])
    
    # sum_{i not in pair} lambda_i tilde_i
    for i in free_indices:
        if i not in tildes_free:
             return None
        rhs_0 -= lambdas[i] * tildes_free[i][0]
        rhs_1 -= lambdas[i] * tildes_free[i][1]
        
    # Matrix M columns are lambda_a and lambda_b.
    # M * [tilde_a^dot; tilde_b^dot] = RHS^dot (vector)
    # M = [ [lambda_a[0], lambda_b[0]], [lambda_a[1], lambda_b[1]] ]
    
    M = matrix(Field, [[lambdas[a][0], lambdas[b][0]], 
                       [lambdas[a][1], lambdas[b][1]]])
                       
    try:
        # Robust check for singularity
        det = M.det()
        if abs(det) < Field(1e-12): 
            return None
            
        sol_0 = M.solve_right(rhs_0)
        sol_1 = M.solve_right(rhs_1)
        
        tildes = {}
        for i in free_indices:
            tildes[i] = tildes_free[i]
            
        tildes[a] = vector(Field, [sol_0[0], sol_1[0]])
        tildes[b] = vector(Field, [sol_0[1], sol_1[1]])
        
        return tildes
    except Exception:
        return None



