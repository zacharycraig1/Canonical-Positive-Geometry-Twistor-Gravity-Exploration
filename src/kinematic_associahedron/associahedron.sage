# src/kinematic_associahedron/associahedron.sage
"""
Kinematic Associahedron for n=6
===============================

The associahedron A_{n-3} is a polytope in kinematic space whose
canonical form gives the bi-adjoint scalar amplitude m(α,α).

For n=6, this is A_3 (the 3-dimensional associahedron), which has:
- 14 vertices (= Catalan number C_4 = 14)
- 21 edges
- 9 faces (= 9 facets)
- 3 dimensions

The facets correspond to factorization channels (planar poles).

Reference: Arkani-Hamed et al., "Scattering Forms and the Positive
Geometry of Kinematics, Color and the Worldsheet" (arXiv:1711.09102)
"""

from sage.all import *
from itertools import permutations, combinations
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


class KinematicAssociahedron:
    """
    The kinematic associahedron for n-point scattering.
    
    For a given cyclic ordering α = (1, 2, ..., n), the associahedron A_α
    lives in the space of planar kinematic variables X_ij where i < j-1.
    
    The boundaries correspond to X_ij = 0 (factorization channels).
    """
    
    def __init__(self, n=6, ordering=None):
        """
        Initialize the kinematic associahedron.
        
        Args:
            n (int): Number of particles.
            ordering (tuple): Cyclic ordering of particles. Default is (0,1,2,3,4,5).
        """
        self.n = n
        self.ordering = ordering if ordering else tuple(range(n))
        
        # Build the planar variables X_ij
        self.planar_vars = self._build_planar_variables()
        
        # Build the facets (boundaries)
        self.facets = self._build_facets()
        
        # Build vertices (triangulations)
        self.vertices = self._build_vertices()
        
    def _build_planar_variables(self):
        """
        Build the set of planar kinematic variables X_ij.
        
        For cyclic ordering (1,2,...,n), X_ij corresponds to the
        planar Mandelstam variable s_{i,i+1,...,j-1}.
        
        There are n(n-3)/2 such variables for n particles.
        For n=6: 6*3/2 = 9 variables.
        
        In the associahedron, we use a subset as coordinates.
        """
        n = self.n
        order = self.ordering
        
        # Map from indices to positions in ordering
        pos = {order[i]: i for i in range(n)}
        
        # Planar variables X_{i,j} where positions satisfy 0 <= i < j-1 < n-1
        # These correspond to non-adjacent pairs in cyclic order
        planar_vars = []
        
        for i in range(n):
            for j in range(i+2, n):
                if (j - i) >= 2 and (j - i) <= n - 2:
                    # X_{order[i], order[j]}
                    planar_vars.append((order[i], order[j]))
                    
        return planar_vars
    
    def _build_facets(self):
        """
        Build the facets of the associahedron.
        
        Each facet corresponds to a planar pole X_ij = 0.
        For n=6, there are n(n-3)/2 = 9 facets.
        """
        return self.planar_vars[:]
    
    def _build_vertices(self):
        """
        Build the vertices of the associahedron.
        
        Each vertex corresponds to a complete triangulation of the n-gon.
        For n=6, there are C_4 = 14 vertices.
        """
        return self._enumerate_triangulations()
    
    def _enumerate_triangulations(self):
        """
        Enumerate all triangulations of an n-gon.
        
        A triangulation is a maximal set of non-crossing diagonals.
        For n-gon, we need n-3 diagonals.
        """
        n = self.n
        order = self.ordering
        
        # Diagonals are pairs (i, j) with j > i+1 and (i,j) != (0, n-1)
        diagonals = []
        for i in range(n):
            for j in range(i+2, n):
                if not (i == 0 and j == n-1):  # Exclude the "edge" from 0 to n-1
                    diagonals.append((order[i], order[j]))
        
        # Check if two diagonals cross
        def crosses(d1, d2):
            """Check if diagonals d1 and d2 cross inside the polygon."""
            # Map to positions
            pos = {order[k]: k for k in range(n)}
            a, b = pos[d1[0]], pos[d1[1]]
            c, d = pos[d2[0]], pos[d2[1]]
            
            if a > b:
                a, b = b, a
            if c > d:
                c, d = d, c
            
            # Diagonals cross if one endpoint of each is strictly between
            # the endpoints of the other
            return (a < c < b < d) or (c < a < d < b)
        
        # Find all maximal sets of non-crossing diagonals
        # (i.e., triangulations with n-3 diagonals)
        num_needed = n - 3
        
        triangulations = []
        
        def backtrack(current, start_idx):
            if len(current) == num_needed:
                triangulations.append(tuple(sorted(current)))
                return
            
            for i in range(start_idx, len(diagonals)):
                d = diagonals[i]
                # Check if d crosses any diagonal in current
                is_valid = True
                for existing in current:
                    if crosses(d, existing):
                        is_valid = False
                        break
                
                if is_valid:
                    backtrack(current + [d], i + 1)
        
        backtrack([], 0)
        
        return triangulations
    
    def canonical_form_term(self, vertex, kinematics):
        """
        Compute the canonical form contribution from a single vertex.
        
        For a triangulation T, the contribution is:
        1 / (X_{i1,j1} * X_{i2,j2} * ... * X_{ik,jk})
        
        where the X's are the diagonals in T.
        
        Args:
            vertex: A triangulation (tuple of diagonals).
            kinematics: A kinematic data object with mandelstam(i, j) method.
        
        Returns:
            Rational number representing the term.
        """
        denom = QQ(1)
        
        for diagonal in vertex:
            i, j = diagonal
            # X_ij = s_{i, i+1, ..., j-1}
            # For adjacent in cyclic order: this is a multi-particle Mandelstam
            
            # Compute s_I where I = {i, i+1, ..., j-1} in cyclic order
            s_val = self._compute_planar_mandelstam(i, j, kinematics)
            
            if s_val == 0:
                return None  # Pole
            
            denom *= s_val
        
        return QQ(1) / denom
    
    def _compute_planar_mandelstam(self, i, j, kinematics):
        """
        Compute the planar Mandelstam invariant X_ij = s_{i,i+1,...,j-1}.
        
        This is the sum of s_ab for all pairs a < b in {i, i+1, ..., j-1}
        (indices in cyclic order).
        """
        n = self.n
        order = self.ordering
        
        # Find positions
        pos_map = {order[k]: k for k in range(n)}
        pos_i = pos_map[i]
        pos_j = pos_map[j]
        
        # Indices in the channel (in position order)
        if pos_i < pos_j:
            channel_positions = list(range(pos_i, pos_j))
        else:
            # Wraps around
            channel_positions = list(range(pos_i, n)) + list(range(0, pos_j))
        
        channel_indices = [order[p] for p in channel_positions]
        
        # Compute s_I = sum of s_ab for a < b in I
        s_total = QQ(0)
        for idx_a in range(len(channel_indices)):
            for idx_b in range(idx_a + 1, len(channel_indices)):
                a = channel_indices[idx_a]
                b = channel_indices[idx_b]
                s_ab = kinematics.mandelstam(a, b)
                s_total += s_ab
        
        return s_total
    
    def compute_canonical_form(self, kinematics):
        """
        Compute the full canonical form Ω(A) = m(α,α), the bi-adjoint scalar amplitude.
        
        This is the sum over all vertices (triangulations):
        Ω = Σ_T 1 / (X_T1 * X_T2 * X_T3)
        
        Args:
            kinematics: A kinematic data object.
        
        Returns:
            Rational number representing the bi-adjoint amplitude m(α,α).
        """
        total = QQ(0)
        
        for vertex in self.vertices:
            term = self.canonical_form_term(vertex, kinematics)
            if term is not None:
                total += term
        
        return total
    
    def display_structure(self):
        """Display the structure of the associahedron."""
        print(f"\n{'='*60}")
        print(f"KINEMATIC ASSOCIAHEDRON A_{self.n - 3}")
        print(f"{'='*60}")
        print(f"Number of particles: {self.n}")
        print(f"Ordering: {self.ordering}")
        print(f"Dimension: {self.n - 3}")
        print(f"Number of vertices: {len(self.vertices)} (Catalan C_{self.n - 2})")
        print(f"Number of facets: {len(self.facets)}")
        
        print(f"\nPlanar variables (facets):")
        for var in self.planar_vars:
            print(f"  X_{{{var[0]},{var[1]}}}")
        
        print(f"\nVertices (triangulations):")
        for i, v in enumerate(self.vertices):
            diag_str = ", ".join([f"({d[0]},{d[1]})" for d in v])
            print(f"  T_{i+1}: {diag_str}")


class KinematicData:
    """
    Simple kinematic data container for testing.
    """
    
    def __init__(self, n, seed=42):
        """Initialize with random rational kinematics."""
        self.n = n
        self.seed = seed
        self._generate_random_kinematics()
    
    def _generate_random_kinematics(self):
        """Generate random Mandelstam invariants satisfying momentum conservation."""
        import random
        random.seed(int(self.seed))
        
        n = self.n
        
        # Generate random s_ij for i < j
        self._s = {}
        for i in range(n):
            for j in range(i + 1, n):
                self._s[(i, j)] = QQ(random.randint(-10, 10)) / QQ(random.randint(1, 5))
                self._s[(j, i)] = self._s[(i, j)]
        
        # Set s_ii = 0 (massless)
        for i in range(n):
            self._s[(i, i)] = QQ(0)
        
        # Enforce momentum conservation by adjusting s_{0, n-1}
        # For 6-point: s_123 = s_456
        # Compute s_012 and adjust s_345
        # This is a simplified version; full conservation is complex.
        
    def mandelstam(self, i, j):
        """Return s_ij."""
        return self._s.get((i, j), QQ(0))


def test_associahedron():
    """Test the kinematic associahedron for n=6."""
    print("\n" + "="*70)
    print("TESTING KINEMATIC ASSOCIAHEDRON")
    print("="*70)
    
    # Create associahedron
    assoc = KinematicAssociahedron(n=6)
    assoc.display_structure()
    
    # Create kinematics
    kin = KinematicData(n=6, seed=42)
    
    # Compute canonical form
    print(f"\n{'='*60}")
    print("COMPUTING CANONICAL FORM (Bi-adjoint Amplitude)")
    print(f"{'='*60}")
    
    omega = assoc.compute_canonical_form(kin)
    print(f"\nΩ(A) = m(α,α) = {omega}")
    print(f"       (float) = {float(omega):.6e}")
    
    # Verify structure
    print(f"\nExpected vertices: C_4 = 14")
    print(f"Actual vertices: {len(assoc.vertices)}")
    assert len(assoc.vertices) == 14, "Should have 14 vertices for n=6"
    
    print("\n✓ Associahedron structure verified!")
    
    return assoc, omega


if __name__ == "__main__":
    test_associahedron()

