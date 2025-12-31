from sage.all import *

def eval_edge_vars_from_spinors(lambdas, tildes, x=None, y=None):
    """
    Computes the edge variables z_{ij} from spinor kinematics.
    
    z_{ij} = ([ij] / <ij>) * <ix><iy> * <jx><jy>
    
    Args:
        lambdas: dictionary of spinor lambda vectors (indexed by 0..n-1)
        tildes: dictionary of spinor lambda_tilde vectors
        x, y: reference spinors (vectors)
        
    Returns:
        Dictionary mapping edge tuples (i,j) or string "z_i_j" to values.
    """
    n = len(lambdas)
    z = {}
    
    # Ensure x, y are vectors
    if x is None or y is None:
        # Default reference spinors? Or raise error?
        # For now, let's assume they are provided.
        raise ValueError("Reference spinors x, y must be provided")
        
    def bracket_angle(i, j):
        # <ij> = det(lambda_i, lambda_j)
        return lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
        
    def bracket_square(i, j):
        # [ij] = det(tilde_i, tilde_j)
        # Note: sign convention might vary.
        return tildes[i][0]*tildes[j][1] - tildes[i][1]*tildes[j][0]
        
    def bracket_ref(i, ref):
        # <i ref>
        return lambdas[i][0]*ref[1] - lambdas[i][1]*ref[0]
    
    for i in range(n):
        for j in range(i + 1, n):
            # z_{ij}
            angle = bracket_angle(i, j)
            square = bracket_square(i, j)
            
            # C_i = <ix><iy>
            Ci = bracket_ref(i, x) * bracket_ref(i, y)
            Cj = bracket_ref(j, x) * bracket_ref(j, y)
            
            # Avoid division by zero if angle is 0 (collinear)
            # In numerical checks, we should pick generic kinematics.
            if angle == 0:
                raise ValueError(f"Singular kinematics: <{i}{j}> = 0")
                
            val = (square / angle) * Ci * Cj
            
            z[(i, j)] = val
            z[f"z_{i}_{j}"] = val
            
    return z
