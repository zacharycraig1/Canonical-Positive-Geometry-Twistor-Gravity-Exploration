from sage.all import *
import sys
import os

# Add src to path
sys.path.append(os.getcwd())

from src.posgeom.toric import compute_lattice_basis, get_toric_exponents
from src.posgeom.forest_polytope import get_forest_exponents

def test_toric_roundtrip():
    print("Testing Toric Roundtrip...")
    
    # 1. Simple Case: Square in 3D
    # (0,0,0), (1,0,0), (0,1,0), (1,1,0) shifted to z=1 plane
    vertices = [
        [0,0,1],
        [1,0,1],
        [0,1,1],
        [1,1,1]
    ]
    
    d, basis, origin = compute_lattice_basis(vertices, saturate=True)
    print(f"Dimension: {d}")
    print(f"Basis:\n{basis}")
    print(f"Origin: {origin}")
    
    exponents = get_toric_exponents(vertices, basis, origin)
    print(f"Exponents: {exponents}")
    
    # Roundtrip check
    for i, u in enumerate(exponents):
        v_recon = origin + vector(ZZ, u) * basis
        v_orig = vector(ZZ, vertices[i])
        if v_recon != v_orig:
            print(f"FAIL: {v_orig} -> {u} -> {v_recon}")
            sys.exit(1)
            
    print("Square test passed.")
    
    # 2. Forest Polytope Case (n=4)
    print("\nTesting n=4 Forest Polytope...")
    roots = [0, 1, 2]
    n = 4
    exponents_n4, _ = get_forest_exponents(n, roots)
    
    # These exponents are the "vertices" of the polytope in edge-space
    d, basis, origin = compute_lattice_basis(exponents_n4, saturate=True)
    print(f"Dimension: {d}")
    
    toric_u = get_toric_exponents(exponents_n4, basis, origin)
    
    for i, u in enumerate(toric_u):
        v_recon = origin + vector(ZZ, u) * basis
        v_orig = vector(ZZ, exponents_n4[i])
        if v_recon != v_orig:
            print(f"FAIL n=4: {v_orig} -> {u} -> {v_recon}")
            sys.exit(1)
            
    print("Forest n=4 test passed.")

if __name__ == "__main__":
    test_toric_roundtrip()






