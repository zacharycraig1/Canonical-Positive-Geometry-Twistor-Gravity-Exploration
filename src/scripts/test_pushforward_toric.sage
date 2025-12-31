import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.canonical_polytope import eval_canonical_form_dual
from src.posgeom.pushforward_toric import toric_canonical_form_square_explicit

def check_toric_square():
    print("Checking Toric/Product Structure (Square)...")
    
    W = vector(QQ, [1, 1, 1])
    
    # Explicit formula from library
    val_explicit = toric_canonical_form_square_explicit(W)
    
    # Polytope Engine
    square_verts = [
        vector(QQ, [0, 0]),
        vector(QQ, [1, 0]),
        vector(QQ, [1, 1]),
        vector(QQ, [0, 1])
    ]
    # Note: Vertex order in list is 00, 10, 11, 01.
    
    # Use explicit triangulation to ensure correctness
    # T1: 0-1-2 (00, 10, 11) -> Lower Right (Wait, 10, 11 is right edge)
    # T2: 0-2-3 (00, 11, 01) -> Upper Left
    # Indices: 0(00), 1(10), 2(11), 3(01).
    
    # Triangulation used in 'toric_canonical_form_square_explicit':
    # S1: (0,0), (1,0), (0,1). Indices: 0, 1, 3.
    # S2: (1,1), (1,0), (0,1). Indices: 2, 1, 3.
    
    tri_matched = [(0, 1, 3), (2, 1, 3)]
    
    val_engine = eval_canonical_form_dual(W, square_verts, triangulation=tri_matched)
    
    print(f"Explicit Formula: {val_explicit}")
    print(f"Polytope Engine:  {val_engine}")
    
    if abs(val_explicit - val_engine) < 1e-10:
        print("SUCCESS: Explicit Square calculation matches Engine.")
    else:
        print("FAILURE: Mismatch.")
        exit(1)

if __name__ == "__main__":
    check_toric_square()
