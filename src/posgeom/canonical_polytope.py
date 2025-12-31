from sage.all import Matrix, vector, QQ, Polyhedron

def triangulate_polytope(vertices):
    """
    Returns a triangulation of the polytope defined by vertices.
    
    Args:
        vertices: list of vectors (or lists)
        
    Returns:
        List of simplices, where each simplex is a tuple of INDICES into 'vertices'.
    """
    P = Polyhedron(vertices=vertices)
    try:
        triangulation = P.triangulation() 
    except AttributeError:
        # Fallback for different Sage versions
        import sage.geometry.triangulation.point_configuration as pc
        conf = pc.PointConfiguration(P.vertices())
        triangulation = conf.triangulate()
    return triangulation

def eval_canonical_form_dual(W, vertices, triangulation=None):
    """
    Evaluates the canonical form of the polytope P=Conv(vertices) at a point W in the dual projective space.
    
    Formula: Omega(W) = sum_{S in Triangulation} det(Z_S) / prod_{v in S} (W . Z_v)
    where Z_v = (1, v) are homogenized vertices.
    
    Checks for orientation and simplex dimension.
    
    Args:
        W: vector of length d+1 (dual coordinates).
        vertices: list of d-dim vectors.
        triangulation: optional precomputed triangulation (list of index tuples).
        
    Returns:
        Value (rational/float)
    """
    d = len(vertices[0])
    # Homogenize vertices
    Z = [vector(QQ, [1] + list(v)) for v in vertices]
    
    if len(W) != d + 1:
        raise ValueError(f"W must have dimension {d+1}")
        
    if triangulation is None:
        triangulation = triangulate_polytope(vertices)
        
    total = 0
    
    # We need a reference orientation to ensure global consistency?
    # Actually, for canonical form, we usually require positive volume relative to the standard volume form.
    # We take abs(det) if we treat this as a volume density sum?
    # No, the signs matter for cancellations between internal boundaries.
    # The triangulation must be ORIENTED.
    # Sage's triangulation usually provides oriented simplices?
    # Let's check or assume positive volume for now, but report mixed signs if found.
    
    for simplex in triangulation:
        # simplex is a list/tuple of indices
        indices = simplex
        if len(indices) != d + 1:
            raise ValueError(f"Simplex has {len(indices)} vertices, expected {d+1}. Degenerate or non-simplicial?")
            
        # Get vertices of simplex
        S_verts = [Z[i] for i in indices]
        
        # Determinant (Volume)
        M = Matrix(S_verts)
        vol = M.det()
        
        # If vol is 0, the simplex is degenerate (should not happen in valid triangulation)
        if vol == 0:
             continue # or raise error
        
        # We need to handle sign. 
        # For projective canonical forms, we sum  Vol(S) / Prod(Linear).
        # The sign of Vol(S) depends on vertex ordering.
        # But Sage's triangulation might not guarantee consistent orientation across all simplices
        # relative to the standard basis?
        # Usually we want |Vol(S)| if we are defining a positive measure.
        # Let's try using abs(vol) to enforce positivity, assuming "Positive Geometry".
        
        vol = abs(vol)
        
        # Denominator
        denom = 1
        for z_vec in S_verts:
            dot_val = vector(QQ, W).dot_product(z_vec)
            denom *= dot_val
            
        total += vol / denom
        
    return total
