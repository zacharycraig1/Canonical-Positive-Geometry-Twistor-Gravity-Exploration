import sys
from sage.all import *

def simplex_canonical_form(W, vertices):
    """
    Computes the canonical form of a simplex defined by 'vertices' at dual point W.
    
    Formula:
      Omega(W) = det(Z) / prod(W . Z_i)
      where Z_i are homogenized vertices.
      
    Args:
        W: Dual vector (dim d+1)
        vertices: List of d vectors (dim d)
        
    Returns:
        Value (rational)
    """
    d = len(vertices[0])
    # Homogenize
    Z = [vector(QQ, [1] + list(v)) for v in vertices]
    
    if len(Z) != d + 1:
        raise ValueError(f"Simplex in d={d} must have {d+1} vertices.")
        
    # Determinant (Volume factor)
    M = Matrix(Z)
    vol = M.det()
    
    # Denominator
    denom = 1
    for z_vec in Z:
        dot_val = vector(QQ, W).dot_product(z_vec)
        denom *= dot_val
        
    if denom == 0:
        raise ValueError("W is on the boundary of the dual simplex (denominator zero).")
        
    return abs(vol) / denom





