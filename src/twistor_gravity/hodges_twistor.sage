#!/usr/bin/env sage
"""
Hodges Formula in Pure Momentum Twistor Variables
==================================================

Implements the Hodges determinant formula for n-point MHV gravity amplitudes
using only momentum twistor variables Z_i ∈ CP^3.

Key formulas:
    M_n^MHV = <12>^8 × det'(Φ) / normalization

where:
    Φ_ij = [ij] / <ij>  (off-diagonal)
    Φ_ii = diagonal element from gauge fixing

The advantage of twistor variables is that:
- Momentum conservation is automatic
- Positivity conditions have geometric meaning in Gr_+(4,n)
- All brackets derive from 4-brackets ⟨ijkl⟩ = det(Z_i, Z_j, Z_k, Z_l)
"""

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))


class HodgesTwistor:
    """
    Computes Hodges MHV gravity amplitude from momentum twistors.
    
    Given n momentum twistors Z_i ∈ QQ^4, computes:
    - Angle brackets <ij> = ε_AB Z_i^A Z_j^B (2-bracket)
    - Square brackets [ij] from 4-brackets via incidence relations
    - Hodges matrix Φ_ij = [ij]/<ij>
    - MHV amplitude M_n = <12>^8 × det'(Φ)
    """
    
    def __init__(self, twistors, negative_helicity=(0, 1)):
        """
        Initialize with momentum twistor data.
        
        Args:
            twistors: List of n 4-vectors Z_i (QQ or RR)
            negative_helicity: Tuple of particle indices with negative helicity
        """
        self.n = len(twistors)
        self.Z = [vector(parent(twistors[0][0]), z) for z in twistors]
        self.neg_hel = negative_helicity
        
        # Precompute brackets
        self._compute_angle_brackets()
        self._compute_four_brackets()
        self._compute_square_brackets()
    
    def _compute_angle_brackets(self):
        """Compute all 2-brackets <ij> = Z_i^0 Z_j^1 - Z_i^1 Z_j^0."""
        n = self.n
        self.angle = {}
        for i in range(n):
            for j in range(n):
                self.angle[(i, j)] = self.Z[i][0] * self.Z[j][1] - self.Z[i][1] * self.Z[j][0]
    
    def _compute_four_brackets(self):
        """Compute all 4-brackets <ijkl> = det(Z_i, Z_j, Z_k, Z_l)."""
        from itertools import combinations
        n = self.n
        self.four_bracket = {}
        
        for indices in combinations(range(n), 4):
            M = matrix([self.Z[k] for k in indices])
            self.four_bracket[indices] = M.det()
    
    def _compute_square_brackets(self):
        """
        Compute [ij] from 4-brackets using incidence relations.
        
        [ij] = <(i-1) i (j-1) j> / (<(i-1) i> <(j-1) j>)
        """
        n = self.n
        self.square = {}
        
        for i in range(n):
            for j in range(n):
                if i == j:
                    self.square[(i, j)] = 0
                    continue
                
                im1 = (i - 1) % n
                jm1 = (j - 1) % n
                
                # 4-bracket <(i-1) i (j-1) j>
                four_br = self.get_four_bracket(im1, i, jm1, j)
                
                # Denominator <(i-1) i> <(j-1) j>
                denom = self.angle[(im1, i)] * self.angle[(jm1, j)]
                
                if denom == 0:
                    self.square[(i, j)] = None
                else:
                    self.square[(i, j)] = four_br / denom
    
    def get_four_bracket(self, i, j, k, l):
        """Get 4-bracket <ijkl> with proper sign from permutation."""
        indices = (i, j, k, l)
        sorted_indices = tuple(sorted(indices))
        
        if sorted_indices not in self.four_bracket:
            return 0
        
        base = self.four_bracket[sorted_indices]
        
        # Compute sign from permutation parity
        inversions = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if indices[a] > indices[b]:
                    inversions += 1
        
        sign = (-1) ** inversions
        return sign * base
    
    def hodges_matrix(self, reference_x=None, reference_y=None):
        """
        Build the n×n Hodges matrix Φ.
        
        Off-diagonal: Φ_ij = [ij] / <ij>  for i ≠ j
        Diagonal: Φ_ii = -Σ_{j≠i} Φ_ij × <j,x><j,y> / (<i,x><i,y>)
        
        where x, y are reference spinors.
        
        Args:
            reference_x: Reference spinor (default: (1, 0))
            reference_y: Reference spinor (default: (0, 1))
        
        Returns:
            n×n Hodges matrix
        """
        n = self.n
        R = parent(self.Z[0][0])
        
        # Default reference spinors
        if reference_x is None:
            reference_x = vector(R, [1, 0])
        if reference_y is None:
            reference_y = vector(R, [0, 1])
        
        Phi = matrix(R, n, n)
        
        # Off-diagonal elements
        for i in range(n):
            for j in range(n):
                if i != j:
                    ang = self.angle[(i, j)]
                    sq = self.square[(i, j)]
                    
                    if ang == 0 or sq is None:
                        # Singular - skip or handle specially
                        Phi[i, j] = 0
                    else:
                        Phi[i, j] = sq / ang
        
        # Diagonal elements (gauge-dependent)
        # Φ_ii = -Σ_{j≠i} Φ_ij × <j,x><j,y> / (<i,x><i,y>)
        for i in range(n):
            lambda_i = vector(R, [self.Z[i][0], self.Z[i][1]])
            ang_ix = lambda_i[0] * reference_x[1] - lambda_i[1] * reference_x[0]
            ang_iy = lambda_i[0] * reference_y[1] - lambda_i[1] * reference_y[0]
            
            if ang_ix == 0 or ang_iy == 0:
                # Reference spinor issue - try different reference
                continue
            
            diag_sum = R(0)
            for j in range(n):
                if j == i:
                    continue
                lambda_j = vector(R, [self.Z[j][0], self.Z[j][1]])
                ang_jx = lambda_j[0] * reference_x[1] - lambda_j[1] * reference_x[0]
                ang_jy = lambda_j[0] * reference_y[1] - lambda_j[1] * reference_y[0]
                
                term = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
                diag_sum -= term
            
            Phi[i, i] = diag_sum
        
        return Phi
    
    def reduced_determinant(self, Phi, deletion_set=(0, 1, 2)):
        """
        Compute reduced determinant of Hodges matrix.
        
        det'(Φ) = det(Φ with rows/cols in deletion_set removed) / normalization
        
        Args:
            Phi: Full Hodges matrix
            deletion_set: Tuple of 3 indices to delete
        
        Returns:
            Reduced determinant
        """
        n = self.n
        deletion_set = sorted(list(deletion_set))
        
        # Extract submatrix
        keep_indices = [i for i in range(n) if i not in deletion_set]
        Phi_reduced = Phi[keep_indices, :][:, keep_indices]
        
        det_reduced = Phi_reduced.det()
        
        # Normalization factor: (<r1 r2><r2 r3><r3 r1>)^2
        r1, r2, r3 = deletion_set
        ang_12 = self.angle[(r1, r2)]
        ang_23 = self.angle[(r2, r3)]
        ang_31 = self.angle[(r3, r1)]
        
        norm_factor = (ang_12 * ang_23 * ang_31)**2
        
        if norm_factor == 0:
            return None
        
        return det_reduced / norm_factor
    
    def mhv_amplitude(self, deletion_set=None):
        """
        Compute the MHV gravity amplitude.
        
        M_n^MHV = <ab>^8 × det'(Φ)
        
        where a, b are the negative helicity particles.
        
        Args:
            deletion_set: Optional tuple of 3 indices to delete (default: (0,1,2))
        
        Returns:
            MHV amplitude value
        """
        if deletion_set is None:
            deletion_set = (0, 1, 2)
        
        # Build Hodges matrix
        Phi = self.hodges_matrix()
        
        # Compute reduced determinant
        det_prime = self.reduced_determinant(Phi, deletion_set)
        
        if det_prime is None:
            return None, "singular"
        
        # Helicity factor <ab>^8
        a, b = self.neg_hel
        helicity_factor = self.angle[(a, b)]**8
        
        return helicity_factor * det_prime, "ok"
    
    def is_positive(self, verbose=False):
        """
        Check if momentum twistors are in positive region.
        
        Positive means: all ordered 4-brackets <i i+1 j j+1> > 0
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
                
                bracket = self.get_four_bracket(i, ip1, j, jp1)
                
                if bracket <= 0:
                    if verbose:
                        print(f"Negative bracket: <{i} {ip1} {j} {jp1}> = {bracket}")
                    all_positive = False
        
        return all_positive


def test_hodges_twistor():
    """Test HodgesTwistor implementation."""
    print("="*60)
    print("TESTING HODGES IN TWISTOR VARIABLES")
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
    
    # Create HodgesTwistor
    ht = HodgesTwistor(Z, negative_helicity=(0, 1))
    
    # Check positivity
    is_pos = ht.is_positive(verbose=True)
    print(f"\nTwistors in positive region: {is_pos}")
    
    # Compute amplitude
    amp, status = ht.mhv_amplitude()
    print(f"\nMHV amplitude: {amp}")
    print(f"Status: {status}")
    
    # Compare with spinor-based Hodges
    print("\n" + "-"*40)
    print("Comparing with spinor-based Hodges...")
    
    try:
        from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
        
        # Convert twistors to spinors
        lambdas = [vector(QQ, [z[0], z[1]]) for z in Z]
        
        # Reconstruct tilde_lambdas from twistors
        tilde_lambdas = []
        for i in range(n):
            im1 = (i - 1) % n
            ip1 = (i + 1) % n
            
            mu_im1 = vector(QQ, [Z[im1][2], Z[im1][3]])
            mu_i = vector(QQ, [Z[i][2], Z[i][3]])
            mu_ip1 = vector(QQ, [Z[ip1][2], Z[ip1][3]])
            
            ang_i_ip1 = ht.angle[(i, ip1)]
            ang_ip1_im1 = ht.angle[(ip1, im1)]
            ang_im1_i = ht.angle[(im1, i)]
            
            denom = ang_im1_i * ang_i_ip1
            if denom == 0:
                print(f"  Singular at particle {i}")
                tilde_lambdas.append(vector(QQ, [0, 0]))
            else:
                num = mu_im1 * ang_i_ip1 + mu_i * ang_ip1_im1 + mu_ip1 * ang_im1_i
                tilde_lambdas.append(num / denom)
        
        amp_spinor, status_spinor = hodges_npt_mhv_canonical(
            lambdas, tilde_lambdas, (0, 1)
        )
        
        print(f"Spinor-based Hodges: {amp_spinor}")
        
        if amp is not None and amp_spinor is not None:
            ratio = amp / amp_spinor
            print(f"Ratio (twistor/spinor): {float(ratio):.6f}")
    except Exception as e:
        print(f"  Could not compare: {e}")
    
    return ht


if __name__ == "__main__":
    test_hodges_twistor()

