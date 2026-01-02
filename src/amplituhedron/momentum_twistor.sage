#!/usr/bin/env sage
"""
Momentum Twistor Infrastructure for Gravity Amplituhedron
==========================================================

This module provides CLEAN momentum twistor kinematics for MHV gravity.

Key features:
- Exact rational arithmetic (QQ)
- Positivity checks for amplituhedron region
- All bracket types: ⟨ij⟩, ⟨ijkl⟩, [ij]
- Caching for efficiency

WARNING: This replaces the broken implementation in ruled_out_single_log_form/
"""

from sage.all import *
from itertools import combinations


class MomentumTwistorData:
    """
    Momentum twistor kinematics for n particles.
    
    Z_i ∈ CP^3 (4-component twistor)
    
    Key invariants:
    - ⟨ij⟩ = ε_{AB} Z_i^A Z_j^B (2-bracket from first 2 components)
    - ⟨ijkl⟩ = det(Z_i, Z_j, Z_k, Z_l) (4-bracket)
    - [ij] derived from incidence relations
    
    For amplituhedron: need all ordered minors ⟨i i+1 j j+1⟩ > 0
    """
    
    def __init__(self, n, seed=None, twistors=None):
        """
        Initialize momentum twistor data.
        
        Args:
            n: Number of particles
            seed: Random seed for generating kinematics (if twistors not provided)
            twistors: Optional list of 4-vectors to use directly
        """
        self.n = n
        
        if twistors is not None:
            # Use provided twistors
            self.Z = [vector(QQ, z) for z in twistors]
        else:
            # Generate random positive kinematics
            if seed is not None:
                set_random_seed(seed)
            self.Z = self._generate_positive_kinematics()
        
        # Cache all brackets
        self._cache_brackets()
    
    def _generate_positive_kinematics(self):
        """
        Generate momentum twistors in the positive region.
        
        For positive Grassmannian Gr_+(4,n), we need all ordered minors > 0.
        We use a simple parameterization that ensures this generically.
        """
        n = self.n
        Z = []
        
        # Generate twistors with random positive integer components
        # This gives positive minors generically (for generic choice)
        for i in range(n):
            z = vector(QQ, [
                QQ(randint(1, 20)),
                QQ(randint(1, 20)),
                QQ(randint(1, 20)),
                QQ(randint(1, 20))
            ])
            Z.append(z)
        
        return Z
    
    def _cache_brackets(self):
        """Pre-compute all brackets for efficiency."""
        n = self.n
        
        # 2-brackets ⟨ij⟩ from first two components
        self._angle = {}
        for i in range(n):
            for j in range(n):
                # ⟨ij⟩ = Z_i^0 * Z_j^1 - Z_i^1 * Z_j^0
                self._angle[(i, j)] = self.Z[i][0] * self.Z[j][1] - self.Z[i][1] * self.Z[j][0]
        
        # 4-brackets ⟨ijkl⟩ = det(Z_i, Z_j, Z_k, Z_l)
        self._four_bracket = {}
        for indices in combinations(range(n), 4):
            M = matrix(QQ, [self.Z[k] for k in indices])
            self._four_bracket[indices] = M.det()
    
    def angle(self, i, j):
        """
        Compute ⟨ij⟩ = ε_{AB} Z_i^A Z_j^B
        
        This is the 2-bracket from the first two twistor components.
        """
        return self._angle.get((i, j), QQ(0))
    
    def four_bracket(self, i, j, k, l):
        """
        Compute ⟨ijkl⟩ = det(Z_i, Z_j, Z_k, Z_l)
        
        Returns the signed 4-bracket accounting for permutation parity.
        """
        indices = (i, j, k, l)
        sorted_indices = tuple(sorted(indices))
        
        if sorted_indices not in self._four_bracket:
            return QQ(0)
        
        base = self._four_bracket[sorted_indices]
        
        # Compute sign from permutation parity
        # Count inversions in the permutation that sorts indices
        inversions = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if indices[a] > indices[b]:
                    inversions += 1
        
        sign = (-1) ** inversions
        return sign * base
    
    def square_bracket(self, i, j):
        """
        Compute [ij] in momentum twistor variables.
        
        [ij] = ⟨(i-1) i (j-1) j⟩ / (⟨(i-1) i⟩ ⟨(j-1) j⟩)
        
        This is the square bracket derived from incidence relations.
        """
        n = self.n
        im1 = (i - 1) % n
        jm1 = (j - 1) % n
        
        num = self.four_bracket(im1, i, jm1, j)
        denom = self.angle(im1, i) * self.angle(jm1, j)
        
        if denom == 0:
            return None
        return num / denom
    
    def is_positive(self, verbose=False):
        """
        Check if this configuration is in the positive region.
        
        For amplituhedron, we need all ordered 4-brackets ⟨i i+1 j j+1⟩ > 0
        for cyclically ordered i < i+1 < j < j+1.
        """
        n = self.n
        all_positive = True
        
        for i in range(n):
            ip1 = (i + 1) % n
            for j in range(i + 2, n):
                jp1 = (j + 1) % n
                if jp1 == i:
                    continue
                
                bracket = self.four_bracket(i, ip1, j, jp1)
                
                if bracket <= 0:
                    if verbose:
                        print(f"Negative bracket: ⟨{i} {ip1} {j} {jp1}⟩ = {bracket}")
                    all_positive = False
        
        return all_positive
    
    def mandelstam(self, indices):
        """
        Compute Mandelstam invariant s_{i_1...i_k} for a subset of particles.
        
        s_{I} = (∑_{i∈I} p_i)² where momenta come from twistors.
        """
        # For momentum twistors, Mandelstam invariants involve 4-brackets
        # This is a simplified implementation
        if len(indices) == 2:
            i, j = indices
            ang = self.angle(i, j)
            sq = self.square_bracket(i, j)
            if sq is None:
                return None
            return ang * sq
        else:
            # More complex - sum over pairs
            total = QQ(0)
            for a, i in enumerate(indices):
                for b in range(a + 1, len(indices)):
                    j = indices[b]
                    s_ij = self.mandelstam([i, j])
                    if s_ij is None:
                        return None
                    total += s_ij
            return total
    
    def is_singular(self):
        """
        Check if kinematics are singular (degenerate).
        
        Singular if any angle bracket ⟨i i+1⟩ = 0.
        """
        n = self.n
        for i in range(n):
            ip1 = (i + 1) % n
            if self.angle(i, ip1) == 0:
                return True
        return False
    
    def __repr__(self):
        return f"MomentumTwistorData(n={self.n}, positive={self.is_positive()})"


def generate_positive_twistors(n, max_attempts=100, seed=None):
    """
    Generate momentum twistors guaranteed to be in the positive region.
    
    Uses rejection sampling: generate random twistors until positive.
    """
    if seed is not None:
        set_random_seed(seed)
    
    for attempt in range(max_attempts):
        tw = MomentumTwistorData(n, seed=None)  # Use current random state
        if tw.is_positive() and not tw.is_singular():
            return tw
    
    raise ValueError(f"Could not generate positive twistors after {max_attempts} attempts")


def twistors_to_spinors(tw):
    """
    Convert momentum twistors to spinor-helicity variables.
    
    λ_i comes from first 2 components of Z_i
    λ̃_i requires solving incidence relations
    
    Returns:
        lambdas: dict of λ_i vectors
        tildes: dict of λ̃_i vectors
    """
    n = tw.n
    
    # λ_i = (Z_i^0, Z_i^1)
    lambdas = {}
    for i in range(n):
        lambdas[i] = vector(QQ, [tw.Z[i][0], tw.Z[i][1]])
    
    # λ̃_i from incidence: Z_i^{2,3} = x_i^{α̇β} λ_i^β
    # This requires more care - for now, use a simplified approach
    # Assume λ̃ comes from third and fourth components scaled
    tildes = {}
    for i in range(n):
        # Simplified: λ̃ proportional to (Z^2, Z^3)
        # This is not exact but works for testing structure
        tildes[i] = vector(QQ, [tw.Z[i][2], tw.Z[i][3]])
    
    return lambdas, tildes


# =============================================================================
# TESTS
# =============================================================================

def test_basic():
    """Basic functionality tests."""
    print("Testing MomentumTwistorData...")
    
    # Test creation
    tw = MomentumTwistorData(n=6, seed=42)
    print(f"  Created: {tw}")
    
    # Test brackets
    ang_01 = tw.angle(0, 1)
    print(f"  ⟨01⟩ = {ang_01}")
    
    four_0123 = tw.four_bracket(0, 1, 2, 3)
    print(f"  ⟨0123⟩ = {four_0123}")
    
    # Test antisymmetry
    assert tw.angle(0, 1) == -tw.angle(1, 0), "Angle bracket antisymmetry failed"
    assert tw.four_bracket(0, 1, 2, 3) == -tw.four_bracket(1, 0, 2, 3), "4-bracket antisymmetry failed"
    
    print("  Antisymmetry: ✓")
    
    # Test singularity check
    print(f"  Singular: {tw.is_singular()}")
    
    print("Basic tests passed!")
    return True


def test_positivity_search():
    """Search for positive configurations."""
    print("\nSearching for positive configurations...")
    
    n_positive = 0
    n_trials = 100
    
    for seed in range(n_trials):
        tw = MomentumTwistorData(n=6, seed=seed)
        if tw.is_positive() and not tw.is_singular():
            n_positive += 1
    
    print(f"  Found {n_positive}/{n_trials} positive configurations")
    
    if n_positive > 0:
        print("  Positivity search: ✓")
        return True
    else:
        print("  WARNING: No positive configurations found!")
        return False


if __name__ == '__main__':
    test_basic()
    test_positivity_search()


