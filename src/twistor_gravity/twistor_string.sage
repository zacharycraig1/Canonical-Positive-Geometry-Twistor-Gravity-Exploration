#!/usr/bin/env sage
"""
Twistor String Formula for MHV Gravity
======================================

Implements Skinner's worldsheet-to-twistor-space formula for gravity amplitudes.

Key Concept:
    In twistor string theory, the worldsheet Σ (genus 0, n marked points) maps
    to twistor space PT = CP^3.
    
    For MHV amplitudes:
    - The map is degree 1 (a line in twistor space)
    - All external twistors Z_i must lie on this line
    - This geometric constraint gives MHV structure

Formula Structure:
    M_n^gravity = ∫_{M_{0,n}} δ-function constraints × integrand
    
    where:
    - δ-functions impose that worldsheet maps to a line in PT
    - Integrand encodes polarization and momentum data

Relation to CHY:
    CHY is the "scattering equation" version where worldsheet localizes
    to (n-3)! points. Twistor string provides complementary perspective.

This is a skeleton implementation for exploration.
"""

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))


class TwistorLine:
    """
    Represents a line in twistor space CP^3.
    
    A line is parameterized by two points Z_a, Z_b:
        L(t) = Z_a + t * Z_b  for t ∈ CP^1
    
    For MHV: all external twistors lie on this line.
    """
    
    def __init__(self, Z_a, Z_b):
        """
        Initialize twistor line from two points.
        
        Args:
            Z_a, Z_b: Two 4-vectors defining the line
        """
        self.Z_a = vector(parent(Z_a[0]), Z_a)
        self.Z_b = vector(parent(Z_b[0]), Z_b)
    
    def point(self, t):
        """Get point on line at parameter t."""
        if t == Infinity:
            return self.Z_b
        return self.Z_a + t * self.Z_b
    
    def contains(self, Z, tolerance=1e-10):
        """
        Check if point Z lies on this line.
        
        Z is on line iff det(Z_a, Z_b, Z, Z') = 0 for any Z'.
        Equivalently, Z is in span(Z_a, Z_b).
        """
        # Check if Z can be written as α*Z_a + β*Z_b
        # This is equivalent to rank([Z_a, Z_b, Z]) = 2
        M = matrix([self.Z_a, self.Z_b, Z])
        
        # For exact arithmetic
        try:
            rank = M.rank()
            return rank <= 2
        except:
            # For numerical
            det3s = []
            for i in range(4):
                for j in range(i+1, 4):
                    for k in range(j+1, 4):
                        M3 = matrix([
                            [self.Z_a[i], self.Z_a[j], self.Z_a[k]],
                            [self.Z_b[i], self.Z_b[j], self.Z_b[k]],
                            [Z[i], Z[j], Z[k]]
                        ])
                        det3s.append(abs(M3.det()))
            return max(det3s) < tolerance


class MHVConstraint:
    """
    Encodes the MHV constraint: all twistors on a line.
    
    For n twistors Z_1, ..., Z_n to lie on a line:
    - Choose any 4 of them: <i j k l> = 0 for all choices
    - Or equivalently: dim(span{Z_i}) = 2
    """
    
    def __init__(self, twistors):
        """
        Initialize with external twistor data.
        
        Args:
            twistors: List of n 4-vectors Z_i
        """
        self.n = len(twistors)
        self.Z = [vector(parent(twistors[0][0]), z) for z in twistors]
    
    def check_mhv_constraint(self, verbose=False):
        """
        Check if all twistors lie on a common line.
        
        Returns:
            (bool, TwistorLine or None)
        """
        # For n twistors to be collinear in PT, rank(Z_1,...,Z_n) = 2
        M = matrix(self.Z)
        
        try:
            rank = M.rank()
        except:
            rank = 4  # Assume generic
        
        if rank > 2:
            if verbose:
                print(f"Rank of twistor matrix: {rank} (need 2 for MHV)")
            return False, None
        
        # Extract the line
        if rank == 2:
            # Find two linearly independent twistors
            Z_a = self.Z[0]
            Z_b = None
            for i in range(1, self.n):
                if matrix([Z_a, self.Z[i]]).rank() == 2:
                    Z_b = self.Z[i]
                    break
            
            if Z_b is not None:
                return True, TwistorLine(Z_a, Z_b)
        
        return rank <= 2, None


class TwistorStringGravity:
    """
    Computes MHV gravity amplitude via twistor string approach.
    
    The key insight is that for MHV:
    - External data defines points on CP^1 (the worldsheet)
    - The map to twistor space is degree 1 (a line)
    - The amplitude is the integral over the moduli of this line
    
    For 6-point MHV gravity:
        M_6 = ∫ d^2 (line data) × δ(scattering eqns) × [integrand]
    """
    
    def __init__(self, twistors, worldsheet_points=None):
        """
        Initialize twistor string computation.
        
        Args:
            twistors: External momentum twistors Z_i
            worldsheet_points: Optional σ_i coordinates on CP^1
        """
        self.n = len(twistors)
        self.Z = [vector(parent(twistors[0][0]), z) for z in twistors]
        
        # Worldsheet coordinates (if provided)
        if worldsheet_points is not None:
            self.sigma = worldsheet_points
        else:
            # Default gauge fixing
            self.sigma = None
        
        # Precompute brackets
        self._compute_brackets()
    
    def _compute_brackets(self):
        """Compute angle brackets."""
        self.angle = {}
        for i in range(self.n):
            for j in range(self.n):
                self.angle[(i, j)] = self.Z[i][0] * self.Z[j][1] - self.Z[i][1] * self.Z[j][0]
    
    def worldsheet_propagator(self, i, j):
        """
        Worldsheet propagator 1/(σ_i - σ_j).
        
        This appears in the MHV vertex rules.
        """
        if self.sigma is None:
            raise ValueError("Worldsheet coordinates not set")
        
        diff = self.sigma[i] - self.sigma[j]
        if diff == 0:
            return None
        return 1 / diff
    
    def mhv_vertex(self, subset_indices):
        """
        Compute MHV vertex contribution for subset of particles.
        
        The MHV vertex is the basic building block in twistor string.
        
        Args:
            subset_indices: Particles attached to this vertex
        
        Returns:
            Vertex contribution
        """
        # For MHV: the vertex is essentially the Parke-Taylor factor
        # <12>^4 / (<12><23>...<n1>)
        
        m = len(subset_indices)
        if m < 3:
            return 1
        
        # Cyclic product of angle brackets
        product = 1
        for k in range(m):
            i = subset_indices[k]
            j = subset_indices[(k + 1) % m]
            bracket = self.angle[(i, j)]
            if bracket == 0:
                return None
            product *= bracket
        
        return 1 / product
    
    def scattering_amplitude(self, method='vertex'):
        """
        Compute MHV gravity amplitude.
        
        Args:
            method: 'vertex' (MHV vertices) or 'direct' (direct computation)
        
        Returns:
            Amplitude value
        """
        if method == 'vertex':
            # Single MHV vertex with all particles
            all_indices = list(range(self.n))
            vertex = self.mhv_vertex(all_indices)
            
            if vertex is None:
                return None
            
            # Helicity factor
            helicity = self.angle[(0, 1)]**8
            
            # For gravity: vertex squared (schematically)
            return helicity * vertex**2
        
        elif method == 'direct':
            # Use Hodges formula (delegate to hodges_twistor)
            from src.twistor_gravity.hodges_twistor import HodgesTwistor
            ht = HodgesTwistor(self.Z)
            return ht.mhv_amplitude()
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def check_line_geometry(self, verbose=False):
        """
        Check the geometric structure of MHV constraint.
        
        For MHV: twistors should span a 2D subspace (a line in PT).
        
        Returns:
            dict with geometric analysis
        """
        M = matrix(self.Z)
        
        try:
            rank = M.rank()
        except:
            rank = min(4, self.n)
        
        result = {
            'twistor_rank': rank,
            'is_mhv': rank == 2,
            'is_generic': rank == min(4, self.n)
        }
        
        if verbose:
            print(f"Twistor configuration:")
            print(f"  Rank: {rank}")
            print(f"  MHV (rank=2): {result['is_mhv']}")
            print(f"  Generic: {result['is_generic']}")
        
        return result


def test_twistor_string():
    """Test TwistorStringGravity implementation."""
    print("="*60)
    print("TESTING TWISTOR STRING GRAVITY")
    print("="*60)
    
    # Generate random momentum twistors
    set_random_seed(42)
    n = 6
    Z = []
    for i in range(n):
        z = vector(QQ, [QQ(randint(1, 10)) for _ in range(4)])
        Z.append(z)
    
    print(f"\nMomentum twistors:")
    for i, z in enumerate(Z):
        print(f"  Z_{i} = {z}")
    
    # Create twistor string
    ts = TwistorStringGravity(Z)
    
    # Check geometry
    print(f"\nGeometric analysis:")
    geometry = ts.check_line_geometry(verbose=True)
    
    # Check MHV constraint
    print(f"\nMHV constraint check:")
    mhv = MHVConstraint(Z)
    is_mhv, line = mhv.check_mhv_constraint(verbose=True)
    print(f"  Satisfies MHV constraint: {is_mhv}")
    
    # Compute amplitude
    print(f"\nAmplitude computation:")
    
    try:
        amp_vertex = ts.scattering_amplitude('vertex')
        print(f"  MHV vertex method: {float(amp_vertex) if amp_vertex else 'singular'}")
    except Exception as e:
        print(f"  MHV vertex method failed: {e}")
    
    try:
        amp_direct = ts.scattering_amplitude('direct')
        print(f"  Direct (Hodges) method: {amp_direct}")
    except Exception as e:
        print(f"  Direct method failed: {e}")
    
    return ts


if __name__ == "__main__":
    test_twistor_string()

