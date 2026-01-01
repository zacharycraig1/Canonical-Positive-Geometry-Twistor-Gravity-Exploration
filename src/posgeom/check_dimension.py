
import sys
import os
sys.path.append(os.getcwd())
from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.intrinsic_lattice import IntrinsicLattice

def check_dim():
    n = 6
    roots = [0, 1, 2]
    exponents, edge_order = get_forest_exponents(n, roots)
    lattice = IntrinsicLattice(exponents)
    print(f"Number of edges: {len(edge_order)}")
    print(f"Lattice rank (affine dimension): {lattice.rank}")
    print(f"Basis size: {lattice.B.nrows()} x {lattice.B.ncols()}")

if __name__ == "__main__":
    check_dim()


