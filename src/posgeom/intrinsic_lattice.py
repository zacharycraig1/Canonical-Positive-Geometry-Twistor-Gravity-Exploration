import sys
import os

# Ensure we can import from src
sys.path.append(os.getcwd())

# Ensure we can import sage
try:
    from sage.all import Matrix, ZZ, vector, RR
except ImportError:
    print("This module must be run with SageMath: 'sage -python ...'")
    sys.exit(1)

class IntrinsicLattice:
    def __init__(self, exponents):
        """
        Computes an integer basis for the lattice spanned by the exponent vectors.
        
        Args:
            exponents: list of integer vectors (list or tuple)
        """
        if not exponents:
            raise ValueError("Exponents list cannot be empty")
            
        self.exponents = [vector(ZZ, v) for v in exponents]
        self.a0 = self.exponents[0]
        
        # Differences from the first exponent (affine anchor)
        diffs = [v - self.a0 for v in self.exponents[1:]]
        
        if not diffs:
            # Single point, dimension 0
            self.B = Matrix(ZZ, len(self.a0), 0)
            self.basis_vectors = []
        else:
            # Create matrix with diffs as rows
            M = Matrix(ZZ, diffs)
            
            # Row Hermite Normal Form
            # The non-zero rows of the HNF form a basis for the lattice spanned by the rows of M.
            # Sage's hermite_form returns H.
            H = M.hermite_form()
            
            # Extract non-zero rows
            self.basis_vectors = []
            for row in H.rows():
                if not row.is_zero():
                    self.basis_vectors.append(row)
            
            # Basis matrix B (columns are basis vectors)
            # B is m x d where m is ambient dimension, d is lattice dimension
            self.B = Matrix(ZZ, self.basis_vectors).transpose()
            
        self.dim = self.B.ncols()
        self.ambient_dim = self.B.nrows()
        
    def to_intrinsic(self, point_ambient):
        """
        Maps a point in the ambient space (minus a0) to intrinsic coordinates t.
        Solves B * t = point_ambient - a0
        """
        diff = vector(point_ambient) - self.a0
        try:
            return self.B.solve_right(diff)
        except ValueError:
            # Not in the column space
            return None
            
    def to_ambient(self, t):
        """
        Maps intrinsic coordinates t to ambient space.
        Returns a0 + B * t
        """
        return self.a0 + self.B * vector(t)
        
    @property
    def covolume(self):
        """
        Returns the covolume of the lattice (volume of fundamental domain).
        For a basis B (m x d), this is sqrt(det(B^T * B)).
        """
        if self.dim == 0:
            return 1
        gram = self.B.transpose() * self.B
        return gram.det().sqrt()

def main():
    import argparse
    from src.posgeom.forest_polytope import get_forest_exponents
    
    parser = argparse.ArgumentParser(description="Compute intrinsic lattice basis for forest polytope")
    parser.add_argument("--n", type=int, default=6, help="Number of points")
    parser.add_argument("--roots", type=int, nargs="+", default=[0, 1, 2], help="Root vertices")
    
    args = parser.parse_args()
    
    print(f"Computing lattice for n={args.n}, roots={args.roots}")
    exponents, _ = get_forest_exponents(args.n, args.roots)
    
    lattice = IntrinsicLattice(exponents)
    
    print(f"Ambient dimension: {lattice.ambient_dim}")
    print(f"Intrinsic dimension (rank): {lattice.dim}")
    print(f"Basis B shape: {lattice.B.nrows()}x{lattice.B.ncols()}")
    print(f"Covolume: {lattice.covolume}")
    
    # Check if a0 is zero (it usually isn't)
    print(f"Anchor a0 (first 5 coords): {lattice.a0[:5]}...")

if __name__ == "__main__":
    main()

