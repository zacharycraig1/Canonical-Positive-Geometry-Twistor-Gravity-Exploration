import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.canonical_polytope import eval_canonical_form_dual, triangulate_polytope

def test_canonical_form_toys():
    print("Testing Canonical Form on Toy Polytopes...")
    
    # 1. 2-Simplex (Triangle) in P^2
    # Vertices: (1,0), (0,1), (0,0) -> Homogenized: (1,1,0), (1,0,1), (1,0,0)
    # Actually canon_polytope expects non-homogenized vertices
    print("\n--- Test 1: 2-Simplex (Triangle) ---")
    triangle_verts = [vector(QQ, [0, 0]), vector(QQ, [1, 0]), vector(QQ, [0, 1])]
    
    # Dual point W (dimension d+1 = 3)
    # W = (W0, W1, W2)
    # Omega = Area / (W.Z0 W.Z1 W.Z2) ? No, standard formula.
    # Vertices Z0=(1,0,0), Z1=(1,1,0), Z2=(1,0,1)
    
    W = vector(QQ, [1, 1, 1])
    
    val1 = eval_canonical_form_dual(W, triangle_verts)
    print(f"Value at W=[1,1,1]: {val1}")
    
    # Check scaling (degree -d-1 = -3)
    W_scale = vector(QQ, [2, 2, 2])
    val_scale = eval_canonical_form_dual(W_scale, triangle_verts)
    ratio = val1 / val_scale
    print(f"Scaling ratio (2x): {ratio} (Expected 8)")
    
    if abs(ratio - 8) > 1e-10:
        print("FAILURE: Scaling property violated.")
        exit(1)
        
    # 2. Square (Triangulation Independence)
    print("\n--- Test 2: Square (Triangulation Independence) ---")
    # Square: (0,0), (1,0), (1,1), (0,1)
    square_verts = [
        vector(QQ, [0, 0]),
        vector(QQ, [1, 0]),
        vector(QQ, [1, 1]),
        vector(QQ, [0, 1])
    ]
    
    # Triangulation A: Split by diagonal (0,0)-(1,1) -> T1: 0-1-2, T2: 0-2-3
    # Indices: 0, 1, 2 and 0, 2, 3
    tri_A = [(0, 1, 2), (0, 2, 3)]
    
    # Triangulation B: Split by diagonal (1,0)-(0,1) -> T1: 0-1-3, T2: 1-2-3
    # Indices: 0, 1, 3 and 1, 2, 3
    tri_B = [(0, 1, 3), (1, 2, 3)]
    
    W_test = vector(QQ, [2, 1, 3])
    
    val_A = eval_canonical_form_dual(W_test, square_verts, triangulation=tri_A)
    val_B = eval_canonical_form_dual(W_test, square_verts, triangulation=tri_B)
    
    print(f"Triangulation A: {val_A}")
    print(f"Triangulation B: {val_B}")
    
    if abs(val_A - val_B) < 1e-10:
        print("SUCCESS: Values match exactly.")
    else:
        print("FAILURE: Triangulation dependence detected!")
        diff = val_A - val_B
        print(f"Difference: {diff}")
        exit(1)
        
    # 3. 3-Simplex (Tetrahedron)
    print("\n--- Test 3: 3-Simplex ---")
    tet_verts = [
        vector(QQ, [0,0,0]),
        vector(QQ, [1,0,0]),
        vector(QQ, [0,1,0]),
        vector(QQ, [0,0,1])
    ]
    W_tet = vector(QQ, [1, 1, 1, 1])
    val_tet = eval_canonical_form_dual(W_tet, tet_verts)
    print(f"Tetrahedron Value: {val_tet}")
    
    print("\nALL TOY TESTS PASSED")

if __name__ == "__main__":
    test_canonical_form_toys()

