from sage.all import *
import sys
import os

sys.path.append(os.getcwd())

from src.posgeom.canonical_polytope import eval_canonical_form_dual
from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.toric import compute_lattice_basis, get_toric_exponents

def test_triangulation_invariance():
    print("Testing Triangulation Invariance (H2)...")
    
    # Use n=4 forest polytope (Square)
    roots = [0, 1, 2]
    n = 4
    ambient_vertices, _ = get_forest_exponents(n, roots)
    d, basis, origin = compute_lattice_basis(ambient_vertices, saturate=True)
    vertices_u = get_toric_exponents(ambient_vertices, basis, origin)
    
    # Vertices of square: (0,0), (1,0), (0,1), (1,1)
    # Order might vary.
    # Let's manually define two triangulations of a square.
    # Square: 0:(0,0), 1:(1,0), 2:(0,1), 3:(1,1)
    # T1: (0,1,3) + (0,3,2)
    # T2: (0,1,2) + (1,3,2)
    
    # Map back to indices in vertices_u list
    # We need to identify indices.
    # Let's force a specific square for testing.
    sq_verts = [[0,0], [1,0], [0,1], [1,1]]
    
    T1 = [(0,1,3), (0,3,2)]
    T2 = [(0,1,2), (1,3,2)]
    
    # Test point W
    W = [1, 2, 3] # Random-ish
    
    val1 = eval_canonical_form_dual(W, sq_verts, triangulation=T1)
    val2 = eval_canonical_form_dual(W, sq_verts, triangulation=T2)
    
    print(f"Val1 (Diag 1): {val1}")
    print(f"Val2 (Diag 2): {val2}")
    
    if abs(val1 - val2) < 1e-9:
        print("SUCCESS: Triangulation invariance holds.")
    else:
        print("FAILURE: Values differ!")
        # Debugging signs
        # T1 S1: (0,0), (1,0), (1,1) -> det((1,0,0),(1,1,0),(1,1,1)) = 1
        # T1 S2: (0,0), (1,1), (0,1) -> det((1,0,0),(1,1,1),(1,0,1)) = -1 -> abs=1
        
        # If we use signed volume, we need consistent orientation.
        # If we use abs volume, we are summing positive contributions?
        # For a CONVEX polytope W inside/outside?
        # Canonical form is rational function.
        # It should be invariant.
        
if __name__ == "__main__":
    test_triangulation_invariance()







