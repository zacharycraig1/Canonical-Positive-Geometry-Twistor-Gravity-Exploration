# src/kinematic_associahedron/klt_geometry.sage
"""
KLT Kernel as Geometric Constraint
==================================

The KLT (Kawai-Lewellen-Tye) relation expresses gravity amplitudes
as a double copy of gauge theory amplitudes:

    M_n^gravity = Σ_{α,β} A_n^YM(α) × S[α|β] × A_n^YM(β)

The KLT kernel S[α|β] is a polynomial in Mandelstam invariants.

Geometrically, the gravity amplitude can be understood as a
canonical form on a "product geometry" weighted by the KLT kernel.

Reference: 
- Kawai, Lewellen, Tye (1986)
- Bern, Carrasco, Johansson (BCJ duality)
- Arkani-Hamed et al. (1711.09102) for geometric interpretation
"""

from sage.all import *
from itertools import permutations
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def klt_kernel_6pt(alpha, beta, kin):
    """
    Compute the field-theory KLT kernel S[α|β] for n=6.
    
    The kernel is defined as:
    S[α|β] = Π_{i∈α'} (s_{1,α_i} + Σ_{j<i, θ_β(α_j,α_i)} s_{α_j,α_i})
    
    where α' = α \ {1, n-1, n} is the permuted subset.
    θ_β(a,b) = 1 if a comes after b in β.
    
    Args:
        alpha: Permutation of {1,2,3} (the permuted subset for n=6)
        beta: Permutation of {1,2,3}
        kin: Kinematic data with mandelstam(i,j) method
    
    Returns:
        Rational value of the kernel
    """
    alpha = list(alpha)
    beta = list(beta)
    
    # Build position map for beta (which element comes after which)
    pos_beta = {beta[i]: i for i in range(len(beta))}
    
    # θ_β(a,b) = 1 if a comes after b in β
    def theta(a, b):
        return 1 if pos_beta.get(a, -1) > pos_beta.get(b, -1) else 0
    
    # Compute the kernel
    kernel = QQ(1)
    
    for i in range(len(alpha)):
        # s_{1, α_i} (where 1 is index 0)
        s_1_ai = kin.mandelstam(0, alpha[i])
        
        term = s_1_ai
        
        # Sum over j < i with θ_β(α_j, α_i) = 1
        for j in range(i):
            if theta(alpha[j], alpha[i]):
                s_aj_ai = kin.mandelstam(alpha[j], alpha[i])
                term += s_aj_ai
        
        kernel *= term
    
    return kernel


def parke_taylor_denominator(order, kin):
    """
    Compute the Parke-Taylor denominator for a given color ordering.
    
    PT(order) = <order[0] order[1]> <order[1] order[2]> ... <order[n-1] order[0]>
    
    For bi-adjoint scalar, this is just the cyclic product of 1's,
    so the amplitude structure comes from the poles only.
    
    Args:
        order: Color ordering as list of indices
        kin: Kinematic data (with angle brackets if needed)
    
    Returns:
        Product of angle brackets, or 1 for scalar theory
    """
    # For scalar theory, the PT factor is trivially 1
    # For YM, we'd use angle brackets
    return QQ(1)


def bi_adjoint_amplitude(alpha_order, beta_order, kin):
    """
    Compute the bi-adjoint scalar amplitude m(α, β).
    
    m(α, β) = Σ_{compatible diagrams} 1 / (X_1 × X_2 × ... × X_{n-3})
    
    The sum is over diagrams (triangulations) compatible with BOTH orderings.
    
    Args:
        alpha_order: First color ordering
        beta_order: Second color ordering
        kin: Kinematic data
    
    Returns:
        Rational amplitude value
    """
    from src.kinematic_associahedron.associahedron import KinematicAssociahedron
    
    n = len(alpha_order)
    
    # Get triangulations for each ordering
    assoc_alpha = KinematicAssociahedron(n, ordering=tuple(alpha_order))
    assoc_beta = KinematicAssociahedron(n, ordering=tuple(beta_order))
    
    # Find compatible triangulations (those in both)
    verts_alpha = set(assoc_alpha.vertices)
    verts_beta = set(assoc_beta.vertices)
    
    common_verts = verts_alpha.intersection(verts_beta)
    
    if len(common_verts) == 0:
        return QQ(0)
    
    # Sum over common triangulations
    total = QQ(0)
    for vertex in common_verts:
        term = assoc_alpha.canonical_form_term(vertex, kin)
        if term is not None:
            total += term
    
    return total


class KLTGeometry:
    """
    The positive geometry for gravity via KLT double copy.
    
    The gravity amplitude is:
    M_n = Σ_{α,β} A(α) × S[α|β] × A(β)
    
    Geometrically, this is related to a "product" of associahedra
    weighted by the KLT kernel.
    """
    
    def __init__(self, n=6):
        """
        Initialize the KLT geometry.
        
        Args:
            n (int): Number of particles.
        """
        self.n = n
        
        # The permuted set (for n=6: {1,2,3} = particles 2,3,4 in 1-indexed)
        self.permuted_set = list(range(1, n-2))  # [1, 2, 3] for n=6
        
        # All permutations of the permuted set
        self.all_perms = list(permutations(self.permuted_set))
        
        # Fixed legs
        self.fixed_first = 0      # Particle 1 (index 0)
        self.fixed_last_2 = n - 2  # Particle n-1 (index 4)
        self.fixed_last = n - 1    # Particle n (index 5)
        
    def compute_klt_kernel_matrix(self, kin):
        """
        Compute the KLT kernel matrix S[α|β] for all permutation pairs.
        
        Returns:
            Matrix of size (n-3)! × (n-3)!
        """
        num_perms = len(self.all_perms)
        S = matrix(QQ, num_perms, num_perms)
        
        for i, alpha in enumerate(self.all_perms):
            for j, beta in enumerate(self.all_perms):
                S[i, j] = klt_kernel_6pt(alpha, beta, kin)
        
        return S
    
    def compute_bi_adjoint_matrix(self, kin):
        """
        Compute the bi-adjoint scalar amplitude matrix m[α|β].
        
        Returns:
            Matrix of size (n-3)! × (n-3)!
        """
        num_perms = len(self.all_perms)
        M = matrix(QQ, num_perms, num_perms)
        
        for i, alpha in enumerate(self.all_perms):
            # Full ordering for α: [0] + list(alpha) + [n-2, n-1]
            alpha_order = [self.fixed_first] + list(alpha) + [self.fixed_last_2, self.fixed_last]
            
            for j, beta in enumerate(self.all_perms):
                # Full ordering for β: [0] + list(beta) + [n-1, n-2] (note swap!)
                beta_order = [self.fixed_first] + list(beta) + [self.fixed_last, self.fixed_last_2]
                
                M[i, j] = bi_adjoint_amplitude(alpha_order, beta_order, kin)
        
        return M
    
    def compute_gravity_amplitude(self, kin, ym_amplitudes_left=None, ym_amplitudes_right=None):
        """
        Compute the gravity amplitude via KLT.
        
        M_grav = Σ_{α,β} A_YM(α) × S[α|β] × A_YM(β)
        
        For pure positive geometry (no helicity), we use m(α,α) as proxy.
        For full MHV gravity, we'd use actual YM amplitudes with helicity.
        
        Args:
            kin: Kinematic data
            ym_amplitudes_left: Vector of A_YM(α) for left orderings (optional)
            ym_amplitudes_right: Vector of A_YM(β) for right orderings (optional)
        
        Returns:
            Gravity amplitude value
        """
        S = self.compute_klt_kernel_matrix(kin)
        
        if ym_amplitudes_left is None or ym_amplitudes_right is None:
            # Use bi-adjoint as proxy (gives scalar gravity-like amplitude)
            M = self.compute_bi_adjoint_matrix(kin)
            
            # For bi-adjoint, the "gravity" is m^T * S * m
            # But actually, for scalar case, it's just tr(M * S)?
            # Let's compute: Σ_{α,β} m(identity, α) * S[α|β] * m(β, identity)
            
            # Use diagonal elements as amplitude vectors
            # A_left = M[identity, :] (row)
            # A_right = M[:, identity] (column)
            
            # Find identity permutation
            identity = tuple(self.permuted_set)
            identity_idx = self.all_perms.index(identity)
            
            A_left = M[identity_idx, :]  # Row vector
            A_right = M[:, identity_idx]  # Column vector
            
            # Gravity = A_left * S * A_right
            return A_left * S * A_right
        else:
            # Use provided YM amplitudes
            A_left = vector(QQ, ym_amplitudes_left)
            A_right = vector(QQ, ym_amplitudes_right)
            
            return A_left * S * A_right
    
    def analyze_geometry(self, kin):
        """
        Analyze the geometric structure of the KLT double copy.
        """
        print(f"\n{'='*70}")
        print("KLT GEOMETRY ANALYSIS")
        print(f"{'='*70}")
        print(f"n = {self.n}")
        print(f"Permuted set: {self.permuted_set}")
        print(f"Number of permutations: {len(self.all_perms)}")
        
        # Compute matrices
        S = self.compute_klt_kernel_matrix(kin)
        M = self.compute_bi_adjoint_matrix(kin)
        
        print(f"\nKLT Kernel Matrix S[α|β]:")
        for i, alpha in enumerate(self.all_perms):
            row_str = " ".join([f"{float(S[i,j]):8.2f}" for j in range(len(self.all_perms))])
            print(f"  {alpha}: [{row_str}]")
        
        print(f"\nBi-adjoint Matrix m[α|β]:")
        for i, alpha in enumerate(self.all_perms):
            row_str = " ".join([f"{float(M[i,j]):8.4f}" for j in range(len(self.all_perms))])
            print(f"  {alpha}: [{row_str}]")
        
        # Check key identities
        print(f"\nKey Identity Check:")
        print(f"  M * S should be related to identity (BCJ duality)")
        
        MS = M * S
        print(f"\n  M * S =")
        for i in range(len(self.all_perms)):
            row_str = " ".join([f"{float(MS[i,j]):8.4f}" for j in range(len(self.all_perms))])
            print(f"    [{row_str}]")
        
        # Compute gravity amplitude
        grav = self.compute_gravity_amplitude(kin)
        print(f"\nGravity amplitude (scalar proxy): {float(grav):.6e}")
        
        return {
            'S_matrix': S,
            'M_matrix': M,
            'gravity_amplitude': grav
        }


def test_klt_geometry():
    """Test the KLT geometry for n=6."""
    from src.kinematic_associahedron.associahedron import KinematicData
    
    print("\n" + "="*70)
    print("TESTING KLT GEOMETRY")
    print("="*70)
    
    # Create geometry
    klt = KLTGeometry(n=6)
    
    # Create kinematics
    kin = KinematicData(n=6, seed=42)
    
    # Analyze
    results = klt.analyze_geometry(kin)
    
    return results


if __name__ == "__main__":
    test_klt_geometry()

