import sys
import os
from sage.all import *
import json

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.forest_polytope_inequalities import get_forest_polytope_inequalities

def compare_facets_n6():
    print("P2: Verifying N=6 Facets against Inequalities...")
    
    n = 6
    roots = [0, 1, 2]
    
    # 1. Get Inequalities
    edges, ineqs = get_forest_polytope_inequalities(n, roots)
    
    # 2. Build Polyhedron in Sage
    # Sage Polyhedron takes ieqs list: [b, -A0, -A1, ...] for b + A.x >= 0 ?
    # Sage convention: [b, a1, a2, ...] corresponds to b + a.x >= 0.
    # Our inequalities: A.x <= b  =>  b - A.x >= 0.
    # So term is b, then -A.
    
    sage_ineqs = []
    for vec, rhs in ineqs:
        # entry = [rhs, -vec[0], -vec[1], ...]
        entry = [rhs] + [-x for x in vec]
        sage_ineqs.append(entry)
        
    print("Constructing Polyhedron (this may take a moment)...")
    P = Polyhedron(ieqs=sage_ineqs)
    
    print(f"Dimension: {P.dim()}")
    print(f"Number of Vertices: {P.n_vertices()}")
    print(f"Number of Facets: {P.n_facets()}")
    
    # Check Dimension
    # Variables: n(n-1)/2 = 15.
    # Equality: 1 (Sum = n-k).
    # Roots R={0,1,2}. Edges between roots must be 0?
    # Yes, x(E(S)) <= |S|-|S cap R|.
    # If S={r1, r2}, |S|=2, |S cap R|=2. x(e) <= 0.
    # So edges between roots are 0.
    # For 3 roots, 3 edges are 0.
    # Expected Dim = 15 - 1 - 3 = 11.
    
    expected_dim = 15 - 1 - 3
    
    if P.dim() == expected_dim:
        print(f"PASS: Dimension is {P.dim()} (Matches rooted forest subspace).")
    else:
        print(f"FAIL: Dimension is {P.dim()} (Expected {expected_dim}).")

    # Facet Analysis
    facets = P.inequalities() # Minimal set
    print(f"Computed Minimal Facets: {len(facets)}")
    
    if P.is_compact():
        print("PASS: Polytope is compact.")
    else:
        print("FAIL: Polytope is unbounded (Error).")
        
    # Check Vertices
    # The number of spanning forests rooted at {0,1,2} in K6.
    # This is det(reduced Laplacian).
    # L_sub = [[6, -1, -1], [-1, 6, -1], [-1, -1, 6]]
    # det = 6(35) - (-7) - (-7)? No.
    # 6*35 + 7 + 7? No.
    # 6*(36-1) - (-1)*(-6-1) + (-1)*(1+6)
    # 210 - 7 - 7 = 196.
    
    # Wait, why did we get 108 vertices?
    # Maybe my manual calculation of det is wrong?
    # Or maybe the polytope vertices are NOT exactly the spanning forests?
    # A matroid base polytope vertices ARE the bases.
    # Let's recompute det with Sage.
    L_sub = Matrix(QQ, [[5, -1, -1], [-1, 5, -1], [-1, -1, 5]])
    expected_vertices = L_sub.det()
    print(f"Expected Vertices (Matrix-Tree): {expected_vertices}")
    
    if P.n_vertices() == expected_vertices:
        print("PASS: Vertex count matches Matrix-Tree Theorem.")
    else:
        print(f"FAIL: Vertex count {P.n_vertices()} != Expected {expected_vertices}.")
        # Debug: 108 = 3 * 36? 
        # eigs of L_sub:
        # J_3 has eigs 3, 0, 0.
        # 6I - J_3 has eigs 6-3, 6-0, 6-0 = 3, 6, 6.
        # det = 3 * 6 * 6 = 108.
        # AHA!
        # My manual calc was wrong.
        # 6(35) - 7 - 7 = 210 - 14 = 196?
        # row 1: 6 * (36-1) = 210.
        # row 2 (col 2): -1 * (-1 * 6 - (-1)*(-1)) = -1 * (-6 - 1) = 7.
        # row 3 (col 3): -1 * ((-1)*(-1) - 6*(-1)) = -1 * (1 + 6) = -7.
        # 210 + 7 - 7 = 210?
        # No.
        # M = [[6, -1, -1], [-1, 6, -1], [-1, -1, 6]]
        # det = 6(35) - (-1)(-7) + (-1)(7)
        # = 210 - 7 - 7 = 196.
        
        # Why is Sage saying 108?
        # Is my L_sub wrong?
        # K6 Laplacian D - A.
        # Deg = 5.
        # L = 5I - J (off diag -1).
        # Not 6I!
        # n=6. Degree is n-1=5.
        # L = [[5, -1, ...], ...]
        # L_sub = [[5, -1, -1], [-1, 5, -1], [-1, -1, 5]]
        # Eigs of J3: 3, 0, 0.
        # Eigs of 5I - J3: 2, 5, 5. (Wait, 5-3=2).
        # det = 2 * 5 * 5 = 50.
        
        # Wait, super-root construction?
        # Vertices of G' are {0..5, rho}.
        # Degree of v in 0..5 is 5 (to others) + 1 (to rho) = 6.
        # Degree of rho is 3 (to 0,1,2).
        # Laplacian of G'.
        # Remove row/col rho.
        # Reduced Laplacian.
        # We want forests rooted at R.
        # This corresponds to spanning trees of G' containing edges (rho, r).
        # This is equivalent to contracting {rho, 0, 1, 2} to a single vertex?
        # No.
        # Matrix-Tree for "All roots connected to supernode".
        # Remove rows/cols for roots?
        # If we remove rows/cols {0, 1, 2}, we get the minor for the remaining vertices {3, 4, 5}.
        # The remaining block is the principal submatrix of L(G') for {3, 4, 5}.
        # Deg(3) = 6? No, 3 is connected to 0,1,2,4,5. Deg=5. Not connected to rho.
        # So diagonal is 5.
        # Off-diagonal is -1.
        # L_sub for {3,4,5} is [[5, -1, -1], [-1, 5, -1], [-1, -1, 5]].
        # det = 50.
        
        # Why did I get 108?
        # Let's check my manual matrix construction logic.
    stats = {
        "n": int(n),
        "roots": [int(r) for r in roots],
        "dim": int(P.dim()),
        "n_vertices": int(P.n_vertices()),
        "n_facets": int(P.n_facets()),
        "expected_dim": int(expected_dim),
        "expected_vertices": int(expected_vertices)
    }
    
    # Save
    with open("RESULTS/facet_certificate_n6.json", "w") as f:
        json.dump(stats, f, indent=2)
        print("Saved certificate.")

if __name__ == "__main__":
    compare_facets_n6()

