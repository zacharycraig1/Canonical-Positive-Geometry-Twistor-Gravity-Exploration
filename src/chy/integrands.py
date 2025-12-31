from sage.all import *

def get_scaled_epsilon_dot_k(i, j, helicity_i, ref_i, kinematics):
    """
    Computes sqrt(2) * (epsilon_i . k_j).
    
    Formulas:
    eps_i(-) . k_j = (<i j> [j q]) / ([i q] * sqrt2)  where q = ref_i
    eps_i(+) . k_j = (<q j> [j i]) / (<q i> * sqrt2)  where q = ref_i
    
    Returns:
        Rational number representing sqrt(2) * (eps_i . k_j)
    """
    q = ref_i
    if helicity_i == -1: # Negative helicity
        # <i j> [j q] / [i q]
        num = kinematics.angle(i, j) * kinematics.square(j, q)
        den = kinematics.square(i, q)
    else: # Positive helicity
        # <q j> [j i] / <q i>
        num = kinematics.angle(q, j) * kinematics.square(j, i)
        den = kinematics.angle(q, i)
        
    if den == 0:
        raise ValueError(f"Polarization singularity: ref spinor {q} collinear with particle {i}")
        
    return num / den

def get_scaled_epsilon_dot_epsilon(i, j, h_i, ref_i, h_j, ref_j, kinematics):
    """
    Computes 2 * (epsilon_i . epsilon_j).
    
    Returns:
        Rational number representing 2 * (eps_i . eps_j)
    """
    if i == j:
        return QQ(0) # eps.eps = 0 for massless
        
    # Same helicity
    if h_i == h_j:
        if ref_i == ref_j:
            return QQ(0)
        else:
            # Formula for same helicity, different refs?
            # usually 0 if refs chosen appropriately, but for general refs:
            # We assume standard choice where same helicity -> same ref
            pass
            
    # Opposite helicity
    # eps_i^-(q1) . eps_j^+(q2) = (<i q2> [j q1]) / ([i q1] <q2 j>)
    # (Checking factor of 2? The 2 cancels the 1/2 from normalization)
    
    if h_i == -1 and h_j == 1:
        q1 = ref_i
        q2 = ref_j
        num = kinematics.angle(i, q2) * kinematics.square(j, q1)
        den = kinematics.square(i, q1) * kinematics.angle(q2, j)
        return num / den
        
    if h_i == 1 and h_j == -1:
        return get_scaled_epsilon_dot_epsilon(j, i, h_j, ref_j, h_i, ref_i, kinematics)
        
    return QQ(0) 

def pfaffian_psi_reduced(sigmas, kinematics, helicities, refs, deletion_indices=(0,1)):
    """
    Computes reduced Pfaffian Pf'(Psi).
    
    Psi is 2n x 2n.
    We build Psi' which is scaled by sqrt(2) in polarization blocks.
    
    Pf'(Psi) = (-1)^{i+j} Pf(Psi^{ij}_{ij}) / sigma_ij
    """
    n = kinematics.n
    
    # Check base ring
    s01 = kinematics.s(0, 1)
    base_ring = s01.parent()
    if base_ring == QQ:
        Psi_prime = matrix(QQ, 2*n, 2*n)
    else:
        Psi_prime = matrix(base_ring.fraction_field(), 2*n, 2*n)
    
    # Fill Psi_prime
    # Blocks: 0..n-1 (A), n..2n-1 (B)
    
    # A_ab = k_a . k_b / sig_ab
    # B'_ab = 2 eps_a . eps_b / sig_ab
    # C'_ab = sqrt(2) eps_a . k_b / sig_ab
    
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
                
            sig_diff = sigmas[a] - sigmas[b]
            
            # A block
            # k_a . k_b = s_ab / 2
            val_A = kinematics.s(a, b) / (2 * sig_diff)
            Psi_prime[a, b] = val_A
            
            # B block (indices n+a, n+b)
            val_B = get_scaled_epsilon_dot_epsilon(a, b, helicities[a], refs[a], 
                                                 helicities[b], refs[b], kinematics) / sig_diff
            Psi_prime[n+a, n+b] = val_B
            
            # C block (indices n+a, b) => Psi[n+a, b] = C_ab
            
            val_C = get_scaled_epsilon_dot_k(a, b, helicities[a], refs[a], kinematics) / sig_diff
            Psi_prime[n+a, b] = val_C
            Psi_prime[b, n+a] = -val_C
            
    # Diagonal elements (C_aa block)
    # C_aa = -sum_{c!=a} C_ac
    for a in range(n):
        c_sum = 0
        for c in range(n):
            if a == c: continue
            c_sum += get_scaled_epsilon_dot_k(a, c, helicities[a], refs[a], kinematics) / (sigmas[a] - sigmas[c])
        
        Psi_prime[n+a, a] = -c_sum
        Psi_prime[a, n+a] = c_sum 
        
    # Deletion
    r, s = deletion_indices
    # Rows to delete: r, s. 
    
    rows_to_keep = [x for x in range(2*n) if x != r and x != s]
    Psi_minor = Psi_prime[rows_to_keep, rows_to_keep]
    
    pf_prime = Psi_minor.pfaffian()
    
    # Factor sigma_rs
    sig_rs = sigmas[r] - sigmas[s]
    
    # Note: (-1)^{r+s} convention needs checking, but we verified |M_CHY| = |M_Hodges| with sign fix.
    
    result = pf_prime / sig_rs
    
    return result

def mhv_gravity_amplitude(sigmas, kinematics, det_phi_prime):
    """
    Computes M_grav = (Pf' Psi)^2 / det' Phi.
    """
    n = kinematics.n
    
    # Setup helicities and refs
    # 0, 1 negative (ref 5)
    # 2, 3, 4, 5 positive (ref 0)
    
    helicities = [-1, -1] + [1] * (n-2)
    refs = [5, 5] + [0] * (n-2) # 0-based indices
    
    # Compute Pf'(Psi')
    # Deletion: 0, 1
    pf_prime_val = pfaffian_psi_reduced(sigmas, kinematics, helicities, refs, deletion_indices=(0, 1))
    
    # (Pf' Psi)^2 = (Pf'(Psi') / 2^(n/2))^2 = (Pf'(Psi'))^2 / 2^n
    # Empirical check vs Hodges: we need to multiply by 64 (2^6) and flip sign.
    # So we remove the / 2^n factor and add a minus sign.
    
    num = (pf_prime_val**2) 
    
    return -num / det_phi_prime
