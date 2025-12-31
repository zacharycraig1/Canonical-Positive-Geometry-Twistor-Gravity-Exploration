from sage.all import *

def detprime_phi(sigmas, kinematics, deletion_set=None):
    """
    Computes reduced determinant det'(Phi) for the CHY matrix Phi.
    
    Phi_ab = s_ab / (sigma_a - sigma_b)^2  for a != b
    Phi_aa = -sum_{b!=a} Phi_ab
    
    det'(Phi) = det(Phi[ijk][ijk]) / (sigma_ij * sigma_jk * sigma_ki)^2
    
    Args:
        sigmas: List of sigma values
        kinematics: SpinorKinematics object (for s_ab)
        deletion_set: Tuple (i, j, k) of indices to delete. Default (0, 1, 2).
        
    Returns:
        det'(Phi) value (QQ)
    """
    n = kinematics.n
    if deletion_set is None:
        deletion_set = (0, 1, 2)
        
    i_idx, j_idx, k_idx = deletion_set
    
    # Build Phi matrix
    # Allow for symbolic entries (Polynomial Ring)
    # Check if we need to infer ring from sigmas/kinematics
    
    # Simple heuristic: try to deduce base ring from input
    # If s_ab is polynomial, use parent of s_ab.
    s01 = kinematics.s(0, 1)
    base_ring = s01.parent()
    if base_ring == QQ: # Default fallback
         Phi = matrix(QQ, n, n)
    else:
         Phi = matrix(base_ring.fraction_field(), n, n)
    
    for a in range(n):
        row_sum = 0
        for b in range(n):
            if a == b:
                continue
                
            s_ab = kinematics.s(a, b)
            diff = sigmas[a] - sigmas[b]
            val = s_ab / (diff**2)
            
            Phi[a, b] = val
            row_sum -= val
            
        Phi[a, a] = row_sum
        
    # Create minor by deleting rows/cols i,j,k
    keep_indices = [x for x in range(n) if x not in deletion_set]
    Phi_minor = Phi[keep_indices, keep_indices]
    
    det0 = Phi_minor.det()
    
    # Compute Vandermonde factor for deletion set
    sig_i = sigmas[i_idx]
    sig_j = sigmas[j_idx]
    sig_k = sigmas[k_idx]
    
    # (sigma_ij * sigma_jk * sigma_ki)^2
    # Note: sigma_ij = sigma_i - sigma_j
    
    vander = ((sig_i - sig_j) * (sig_j - sig_k) * (sig_k - sig_i))**2
    
    if vander == 0:
         raise ValueError("Collision in sigma values for deletion set")
         
    return det0 / vander

def pfaffian_psi_reduced(sigmas, kinematics, pol_vectors, deletion_indices=(0,1)):
    """
    Computes reduced Pfaffian Pf'(Psi).
    Pf'(Psi) = (-1)^{i+j} * Pf(Psi^{ij}_{ij}) / sigma_ij
    
    (Using standard definition from CHY literature, check sign conventions)
    
    Args:
        sigmas: List of sigma values
        kinematics: SpinorKinematics object (k_a)
        pol_vectors: List of polarization vectors (eps_a) (4-vectors or similar)
        deletion_indices: (i, j) to remove
        
    Returns:
        Reduced Pfaffian
    """
    # Placeholder for full implementation if needed.
    # For MHV gravity specialized, we might compute directly from 
    # (Pf'Psi)^2 = det' * M_YM^2 or similar, but primary goal requires exact match.
    # Let's start with det' Phi first as per plan A2.
    pass

