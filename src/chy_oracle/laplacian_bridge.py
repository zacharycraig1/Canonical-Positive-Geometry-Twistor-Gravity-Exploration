import sys
import os
from sage.all import *

try:
    from src.chy_oracle.amplitude_spinor import ang_bracket
except ImportError:
    from chy_oracle.amplitude_spinor import ang_bracket
try:
    from src.chy_oracle.matrix_tree import hodges_weighted_laplacian
except ImportError:
    from chy_oracle.matrix_tree import hodges_weighted_laplacian

def reconstruct_mhv_from_laplacian(lambdas, tilde_lambdas, x, y, roots=(0,1,2), neg=(0,1)):
    """
    Reconstructs the n-point MHV amplitude from the Weighted Laplacian.
    
    Formula:
      M_n = (-1)^(n-1) * <ab>^8 * det(L^(roots)) / ( prod_{k not in roots} C_k^2 * norm(roots) )
    
    where:
      - roots: Indices of deleted rows/cols (size 3)
      - neg: Indices of negative helicity legs (size 2) -> <ab>^8 factor
      - norm(roots): (<r1 r2><r2 r3><r3 r1>)^2
    """
    n = len(lambdas)
    
    if len(roots) != 3:
        raise ValueError("Roots must be a set of size 3 (corank 3)")
        
    # 1. Build Weighted Laplacian
    try:
        L_tilde, C, _ = hodges_weighted_laplacian(lambdas, tilde_lambdas, x, y)
    except ValueError as e:
        return None, str(e)
        
    # 2. Compute Reduced Determinant (Minor)
    # Delete rows/cols in 'roots', keep the rest
    indices_keep = [i for i in range(n) if i not in roots]
    
    det_minor = L_tilde.matrix_from_rows_and_columns(indices_keep, indices_keep).det()
    
    # 3. Product of C_k^2 for kept indices
    prod_C_sq = 1
    for k in indices_keep:
        prod_C_sq *= C[k]**2
        
    if prod_C_sq == 0:
        return None, "domain_violation_C_product_zero"
        
    # 4. Normalization factor for roots
    r1, r2, r3 = roots
    ang_12 = ang_bracket(lambdas[r1], lambdas[r2])
    ang_23 = ang_bracket(lambdas[r2], lambdas[r3])
    ang_31 = ang_bracket(lambdas[r3], lambdas[r1])
    
    norm_factor = (ang_12 * ang_23 * ang_31)**2
    
    if norm_factor == 0:
        return None, "domain_violation_norm_factor_zero"
        
    # 5. Helicity factor
    a, b = neg
    h_factor = ang_bracket(lambdas[a], lambdas[b])**8
    
    # 6. Combine
    # Sign: (-1)^(n-1)
    # For n=6, sign is -1.
    sign = (-1)**(n-1)
    
    M_rec = sign * h_factor * det_minor / (prod_C_sq * norm_factor)
    
    return M_rec, "ok"

def weighted_laplacian(n, lambdas, tildes, x, y):
    """Wrapper for consistency with instructions."""
    return hodges_weighted_laplacian(lambdas, tildes, x, y)




