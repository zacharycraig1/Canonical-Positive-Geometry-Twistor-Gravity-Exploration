import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.pushforward_simplex import simplex_canonical_form
from src.posgeom.canonical_polytope import eval_canonical_form_dual

def check_simplex_pushforward():
    print("Checking Simplex Pushforward Formula...")
    
    # 2-Simplex (Triangle)
    verts = [vector(QQ, [0, 0]), vector(QQ, [1, 0]), vector(QQ, [0, 1])]
    W = vector(QQ, [1, 1, 1])
    
    val_formula = simplex_canonical_form(W, verts)
    val_poly = eval_canonical_form_dual(W, verts)
    
    print(f"Formula: {val_formula}")
    print(f"Polytope: {val_poly}")
    
    if abs(val_formula - val_poly) < 1e-10:
        print("SUCCESS: Simplex formula matches polytope engine.")
    else:
        print("FAILURE: Mismatch.")
        exit(1)
        
    # 3-Simplex (Tetrahedron)
    tet_verts = [
        vector(QQ, [0,0,0]),
        vector(QQ, [1,0,0]),
        vector(QQ, [0,1,0]),
        vector(QQ, [0,0,1])
    ]
    W_tet = vector(QQ, [1, 2, 3, 4])
    
    val_f = simplex_canonical_form(W_tet, tet_verts)
    val_p = eval_canonical_form_dual(W_tet, tet_verts)
    
    print(f"Tetrahedron Formula: {val_f}")
    print(f"Tetrahedron Polytope: {val_p}")
    
    if abs(val_f - val_p) < 1e-10:
        print("SUCCESS: Tetrahedron matches.")
    else:
        print("FAILURE: Mismatch.")
        exit(1)

if __name__ == "__main__":
    check_simplex_pushforward()




