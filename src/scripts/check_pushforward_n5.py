from sage.all import vector, ZZ
import sys
import os

# Add src to path
sys.path.append(os.getcwd())

from src.posgeom.canonical_polytope import eval_canonical_form_dual
from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.toric import compute_lattice_basis, get_toric_exponents

def check_pushforward_n5():
    """
    Check A.4/A.6 for n=5.
    """
    print("Checking n=5 Pushforward (Toric Geometry)...")
    sys.stdout.flush()
    
    # 1. Get Vertices (Exponents)
    roots = [0, 1, 2]
    n = 5
    ambient_vertices, _ = get_forest_exponents(n, roots)
    
    # 2. Project to Affine Basis (Toric Coordinates)
    d, basis, origin = compute_lattice_basis(ambient_vertices, saturate=True)
    print(f"Dimension: {d}")
    
    vertices_u = get_toric_exponents(ambient_vertices, basis, origin)
    print(f"Number of Vertices: {len(vertices_u)}")
    # print(f"Toric Vertices (u): {vertices_u}") # Too many?
    
    # 3. Pick random W in dual space P^{d}
    # Avoid poles! W.Z != 0 for all vertices Z.
    # Vertices have coords 0, 1 mostly.
    # W should be weird.
    import random
    random.seed(42)
    W_test = [random.randint(1, 100) for _ in range(d+1)]
    # Ensure no zeros in scalar products
    
    # This involves triangulating a 6D (?) polytope.
    # Sage might be slow.
    print("Evaluating Canonical Form (this implies triangulation)...")
    val = eval_canonical_form_dual(W_test, vertices_u)
    print(f"Omega({W_test}) = {val}")
    
    if val == 0:
        print("WARNING: Canonical form value is 0!")
    else:
        print("Success: Canonical form evaluated to non-zero value.")

if __name__ == "__main__":
    check_pushforward_n5()

