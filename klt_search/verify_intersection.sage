
from sage.all import *
import sys
import os
sys.path.append(os.getcwd())

load("src/hodges.sage")
from itertools import permutations

# =============================================================================
# COHOMOLOGY BASIS
# =============================================================================

def get_basis_forms(n):
    """
    Generate the (n-3)! Parke-Taylor forms for the basis.
    Basis: permutations of {2, ..., n-2} (indices 1 ... n-3).
    Standard choice: Fixed 1, n-1, n. Permute 2...n-2.
    For n=6: Fixed 1, 5, 6. Permute 2, 3, 4.
    """
    # Indices 1, 2, 3 correspond to particles 2, 3, 4
    permuted = [1, 2, 3]
    return sorted(list(permutations(permuted)))

# =============================================================================
# KLT KERNEL (Standard Pivot 1)
# =============================================================================

def klt_kernel_S3(alpha, beta, twistor):
    # Standard field theory limit kernel
    # S[alpha|beta]
    # alpha, beta are perms of {1, 2, 3} (indices)
    
    pos_beta = {val: idx for idx, val in enumerate(beta)}
    def theta(a, b): return 1 if pos_beta[a] > pos_beta[b] else 0
    def s(i, j): return mandelstam_invariant(twistor, i, j)
    
    i1, i2, i3 = alpha
    
    # Formula for S[i1 i2 i3 | j1 j2 j3] with pivot 0 (particle 1)
    # factor 1: s_{0, i1} + theta(i1, i2)s_{i1, i2} + theta(i1, i3)s_{i1, i3}
    t1 = s(0, i1)
    if t1 is None: return 0
    
    if theta(i1, i2): 
        val = s(i1, i2)
        if val is None: return 0
        t1 += val
    if theta(i1, i3): 
        val = s(i1, i3)
        if val is None: return 0
        t1 += val
    
    # factor 2: s_{0, i2} + theta(i2, i3)s_{i2, i3}
    t2 = s(0, i2)
    if t2 is None: return 0
    if theta(i2, i3): 
        val = s(i2, i3)
        if val is None: return 0
        t2 += val
    
    # factor 3: s_{0, i3}
    t3 = s(0, i3)
    if t3 is None: return 0
    
    return t1 * t2 * t3

# =============================================================================
# INTERSECTION MATRIX (Bi-adjoint Scalar)
# =============================================================================

def biadjoint_scalar_amplitude(alpha, beta, twistor):
    """
    Compute m[alpha|beta] = sum_{all planar cubic graphs} (1/prod s_I) * N(alpha) * N(beta)
    Actually, simpler: m = S^{-1}.
    But we want to VERIFY m * S = I.
    
    We need an independent way to compute m.
    m[alpha|beta] is the bi-adjoint scalar amplitude m(1, alpha, n-1, n | 1, beta, n-1, n).
    
    For n=6, this is sum over 14 poles.
    Or we can use the Feynman diagram sum.
    
    Or use the recursive Berends-Giele relation for m.
    """
    # Placeholder: Computing m is hard without a full library.
    # We will compute S, invert it, and check properties of m = S^-1.
    return None

def verify_intersection_property():
    print("Verifying Intersection Matrix Properties (S^-1)...")
    
    # Generate generic point
    # Use integer twistors for exactness
    twistor = MomentumTwistor(n=6, seed=42)
    # Check domain
    if not twistor.domain_ok:
        print("Twistor domain check failed - trying new seed")
        twistor = MomentumTwistor(n=6, seed=43)
        
    print(f"Domain ok: {twistor.domain_ok}")
    
    # Compute S matrix
    basis = get_basis_forms(6)
    dim = len(basis)
    S = matrix(QQ, dim, dim)
    
    print(f"Basis size: {dim}x{dim}")
    
    for r, alpha in enumerate(basis):
        for c, beta in enumerate(basis):
            S[r, c] = klt_kernel_S3(list(alpha), list(beta), twistor)
            
    if S.rank() != dim:
        print(f"[FAIL] S is rank {S.rank()}, expected {dim}")
        return
        
    print("[SUCCESS] S is full rank.")
    
    # Compute m = S^-1
    m = S.inverse()
    
    # Check if m is symmetric (it should be for the canonical intersection form)
    is_sym = m.is_symmetric()
    print(f"Is m = S^-1 symmetric? {is_sym}")
    
    if not is_sym:
        # Check S symmetry? S is generally NOT symmetric.
        # But m should be symmetric if it's an intersection matrix.
        print("S symmetry:", S.is_symmetric())
        print("Deviation from symmetry (m - m.T norm):")
        diff = m - m.transpose()
        # print(diff)
        
    # Check scaling of m
    # S scales as s^3 (mass dim 6)
    # m scales as s^-3 (mass dim -6) -> correct for amplitude
    
    print("\nIntersection Matrix m (first 3x3 block):")
    print(m[:3, :3])

if __name__ == "__main__":
    verify_intersection_property()

