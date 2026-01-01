
import sys
import os
from sage.all import *

# Add project root to path
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor, MomentumTwistor

def sample_positive_twistor(n, seed=None):
    """
    Generates a MomentumTwistor in the positive region using the Moment Curve.
    Z_i = (1, t_i, t_i^2, t_i^3) with 0 < t_1 < ... < t_n.
    We add a small random perturbation to avoid high symmetry if needed,
    but keeping ordered minors positive.
    Actually, Moment Curve is sufficient for positivity.
    To be safe, we pick random t_i in increasing order.
    """
    if seed is not None:
        set_random_seed(seed)
        
    # Pick n sorted random numbers in (0, 10)
    # Using larger range helps separation
    ts = sorted([QQ.random_element(0, 100)/10 for _ in range(n)])
    
    # Ensure distinct
    for i in range(1, n):
        if ts[i] <= ts[i-1]:
            ts[i] = ts[i-1] + QQ(1)/10
            
    Z = []
    for t in ts:
        # Moment curve
        vec = vector(QQ, [1, t, t**2, t**3])
        # Randomize basis? 
        # The positive Grassmannian is invariant under GL(4)+ action?
        # Z -> Z M.
        # But for now, standard frame is fine.
        Z.append(vec)
        
    return MomentumTwistor(n=n, seed=seed, Z=Z)

def probe_limit(n, limit_type, epsilon, positive=False, **kwargs):
    """
    Generates on-shell kinematics near a specific singularity.
    positive (bool): If True, start from positive kinematics.
    """
    if limit_type == 'soft':
        return probe_soft_limit(n, epsilon, positive=positive, **kwargs)
    elif limit_type == 'collinear':
        return probe_collinear_limit(n, epsilon, positive=positive, **kwargs)
    elif limit_type == 'factorization':
        return probe_factorization_limit(n, epsilon, positive=positive, **kwargs)
    else:
        raise ValueError(f"Unknown limit type: {limit_type}")

def probe_soft_limit(n, epsilon, s_idx=None, seed=None, positive=False):
    """
    Soft limit for particle s_idx.
    """
    if s_idx is None:
        s_idx = n - 1
    
    # Generate hard kinematics for n-1 particles
    n_hard = n - 1
    if positive:
        # Use positive twistor for hard part
        twistor_hard = sample_positive_twistor(n_hard, seed=seed)
        lambdas_hard = [twistor_hard.get_lambda(i) for i in range(n_hard)]
        tildes_hard = []
        for i in range(n_hard):
            tildes_hard.append(twistor_hard.get_tilde_lambda(i))
    else:
        lambdas_hard, tildes_hard = sample_spinors_from_twistor(seed=seed, n=n_hard)
    
    if seed is not None:
        set_random_seed(seed)
        
    s_lambda = vector(QQ, [QQ.random_element(), QQ.random_element()])
    while s_lambda == 0:
        s_lambda = vector(QQ, [QQ.random_element(), QQ.random_element()])
        
    s_tilde = vector(QQ, [QQ.random_element(), QQ.random_element()])
    while s_tilde == 0:
        s_tilde = vector(QQ, [QQ.random_element(), QQ.random_element()])
    
    # Recoil logic
    ts = s_tilde
    t0 = tildes_hard[0]
    t1 = tildes_hard[1]
    
    det_01 = t0[0]*t1[1] - t0[1]*t1[0]
    if det_01 == 0:
        raise ValueError("Degenerate hard kinematics [0 1] = 0")
        
    bracket_s1 = ts[0]*t1[1] - ts[1]*t1[0]
    bracket_s0 = ts[0]*t0[1] - ts[1]*t0[0]
    
    a = -epsilon * bracket_s1 / det_01
    b = epsilon * bracket_s0 / det_01
    
    final_lambdas = []
    final_tildes = []
    
    hard_idx = 0
    for i in range(n):
        if i == s_idx:
            final_lambdas.append(epsilon * s_lambda)
            final_tildes.append(s_tilde)
        else:
            l = lambdas_hard[hard_idx]
            t = tildes_hard[hard_idx]
            
            if hard_idx == 0:
                l = l + a * s_lambda
            elif hard_idx == 1:
                l = l + b * s_lambda
                
            final_lambdas.append(l)
            final_tildes.append(t)
            hard_idx += 1
            
    return SpinorKinematics(n, final_lambdas, final_tildes)

def probe_collinear_limit(n, epsilon, i=None, j=None, seed=None, positive=False):
    """
    Collinear limit i || j.
    """
    if i is None: i = n-2
    if j is None: j = n-1
    if seed is not None:
        set_random_seed(seed)
        
    if positive:
        twistor_full = sample_positive_twistor(n, seed=seed)
    else:
        twistor_full = MomentumTwistor(n=n, seed=seed)
    
    # Modify Z_j to be Z_i + epsilon * W
    # If positive, we should ensure the perturbation W keeps us "positive" direction?
    # Or just generic perturbation is fine, as long as we approach the boundary.
    # For Active Facet check, generic approach is better.
    
    # However, if we start deep in positive region, we approach from inside.
    
    W = vector(QQ, [QQ.random_element(), QQ.random_element(), QQ.random_element(), QQ.random_element()])
    
    new_Z = list(twistor_full.Z)
    new_Z[j] = new_Z[i] + epsilon * W
    
    return sample_spinors_from_twistor_custom(new_Z, n)

def probe_factorization_limit(n, epsilon, seed=None, positive=False):
    return None

def sample_spinors_from_twistor_custom(Z_list, n):
    """
    Helper to create spinor kinematics from explicit Z list.
    """
    lambdas = [vector(QQ, [z[0], z[1]]) for z in Z_list]
    
    def get_angle(i, j):
        return Z_list[i][0]*Z_list[j][1] - Z_list[i][1]*Z_list[j][0]
        
    tilde_lambdas = []
    for i in range(n):
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        
        mu_i = vector(QQ, [Z_list[i][2], Z_list[i][3]])
        mu_im1 = vector(QQ, [Z_list[im1][2], Z_list[im1][3]])
        mu_ip1 = vector(QQ, [Z_list[ip1][2], Z_list[ip1][3]])
        
        ang_i_ip1 = get_angle(i, ip1)
        ang_ip1_im1 = get_angle(ip1, im1)
        ang_im1_i = get_angle(im1, i)
        
        denom = ang_im1_i * ang_i_ip1
        if denom == 0:
            raise ValueError(f"Singularity in reconstruction at index {i}")
            
        num = mu_im1 * ang_i_ip1 + mu_i * ang_ip1_im1 + mu_ip1 * ang_im1_i
        tilde_lambdas.append(num / denom)
        
    return SpinorKinematics(n, lambdas, tilde_lambdas)
