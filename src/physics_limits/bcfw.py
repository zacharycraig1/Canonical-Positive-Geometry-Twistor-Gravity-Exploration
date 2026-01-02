import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor

def bcfw_shift_spinors(lambdas, tilde_lambdas, a, b, z):
    """
    Apply BCFW shift on legs a, b:
    lambda_a_hat = lambda_a + z * lambda_b
    tilde_lambda_b_hat = tilde_lambda_b - z * tilde_lambda_a
    
    All other spinors remain unchanged.
    
    Args:
        lambdas: List of lambda spinors (vectors)
        tilde_lambdas: List of tilde_lambda spinors (vectors)
        a, b: Indices of the shifted legs (0-based)
        z: The complex shift parameter
        
    Returns:
        (shifted_lambdas, shifted_tilde_lambdas)
    """
    n = len(lambdas)
    new_lambdas = [v for v in lambdas]
    new_tilde_lambdas = [v for v in tilde_lambdas]
    
    # \hat\lambda_a = \lambda_a + z \lambda_b
    new_lambdas[a] = lambdas[a] + z * lambdas[b]
    
    # \hat\tilde\lambda_b = \tilde\lambda_b - z \tilde\lambda_a
    new_tilde_lambdas[b] = tilde_lambdas[b] - z * tilde_lambdas[a]
    
    return new_lambdas, new_tilde_lambdas

def get_channel_s(lambdas, tilde_lambdas, indices):
    """
    Compute the Mandelstam invariant s_I = (Sum_{i in I} p_i)^2
    
    Args:
        lambdas, tilde_lambdas: Spinors
        indices: List of indices for the channel
        
    Returns:
        s_I value
    """
    # Calculate total momentum P = sum |i>[i|
    P_matrix = Matrix(QQ, 2, 2, 0)
    for i in indices:
        # |i>[i| is a 2x2 matrix
        # lambda = (l0, l1)^T
        # tilde_lambda = (lt0, lt1)
        # |i>[i| = [[l0*lt0, l0*lt1], [l1*lt0, l1*lt1]]
        
        l = lambdas[i]
        lt = tilde_lambdas[i]
        
        term = Matrix(QQ, 2, 2, [
            [l[0]*lt[0], l[0]*lt[1]],
            [l[1]*lt[0], l[1]*lt[1]]
        ])
        P_matrix += term
        
    # P^2 = det(P) (up to factor of 2 or sign depending on metric conventions, usually det(P) = p^2 for SL(2,C))
    # For p_{alpha dot_alpha}, p^2 = det(p).
    # Proof: p_{ab} = v_mu sigma^mu_ab. det(p) = v^2.
    return P_matrix.det()

def solve_bcfw_pole(lambdas, tilde_lambdas, a, b, channel_indices):
    """
    Solve for z_star such that s_channel(z_star) = 0.
    
    Args:
        lambdas, tilde_lambdas: Initial spinors
        a, b: Shift indices
        channel_indices: Indices defining the channel s_I
        
    Returns:
        z_star (or None if no solution/constant)
    """
    # Check linearity:
    # s(z) = s(0) + z * (s(1) - s(0))
    # This is true for BCFW shift on (a, b) for any channel I.
    
    def eval_s(z_val):
        L, Lt = bcfw_shift_spinors(lambdas, tilde_lambdas, a, b, z_val)
        return get_channel_s(L, Lt, channel_indices)
        
    s0 = eval_s(QQ(0))
    s1 = eval_s(QQ(1))
    
    slope = s1 - s0
    
    if slope == 0:
        return None
        
    z_star = -s0 / slope
    
    # Verify
    s_star = eval_s(z_star)
    if abs(s_star) > 1e-10:
        print(f"Warning: solve_bcfw_pole found z_star={z_star} but s(z_star)={s_star}")
        
    return z_star

def get_momentum_matrix(lambdas, tilde_lambdas, indices):
    """
    Compute P_matrix = sum_{i in indices} |i>[i|
    """
    P_matrix = Matrix(QQ, 2, 2, 0)
    for i in indices:
        l = lambdas[i]
        lt = tilde_lambdas[i]
        term = Matrix(QQ, 2, 2, [
            [l[0]*lt[0], l[0]*lt[1]],
            [l[1]*lt[0], l[1]*lt[1]]
        ])
        P_matrix += term
    return P_matrix

def decompose_momentum_spinors(P_matrix, ref_lambda=None):
    """
    Decompose null vector P into lambda, tilde_lambda.
    P = lambda * tilde_lambda
    
    If ref_lambda is provided, we set:
    lambda = ref_lambda
    tilde_lambda = P * ref_lambda / <??> ? No.
    P |mu] = |lambda> [tilde mu].
    
    Wait. The little group scaling is P = (t lambda) * (1/t tilde_lambda).
    
    For BCFW, we usually define lambda_P in terms of the shifted spinors.
    e.g. if we shift a, b, and P is sum of a subset including a (shifted).
    P(z*) = p_a(z*) + ...
    P(z*)^2 = 0.
    
    Standard BCFW Internal Spinor Definition (e.g. Elvang/Huang 2013):
    If P = p_a(z) + K (a is shifted).
    lambda_P = lambda_a(z)
    tilde_lambda_P = P * lambda_ref / <lambda_P lambda_ref> ??
    Wait. P_{ab} = lambda_P_a tilde_lambda_P_b.
    P |xi] = |lambda_P> [tilde_lambda_P xi] ?
    P |mu] = |lambda_P> [tilde_lambda_P mu].
    
    If we pick lambda_P = P |eta], then tilde_lambda_P is fixed.
    
    Args:
        P_matrix: 2x2 matrix P_{alpha, dot_alpha}
        ref_lambda: If not None, we force lambda_P = ref_lambda (must be consistent).
        
    Returns:
        (lambda_vec, tilde_lambda_vec)
    """
    if abs(P_matrix.det()) > 1e-8:
        print(f"Warning: Decomposing non-null momentum, det={P_matrix.det()}")
    
    # If ref_lambda is given, use it.
    if ref_lambda is not None:
        lam = ref_lambda
        # Solve for lt.
        # P = lam * lt. P_ab = lam_a * lt_b.
        # lt_b = P_ab / lam_a.
        if abs(lam[0]) > 1e-10:
            lt = vector([P_matrix[0,0]/lam[0], P_matrix[0,1]/lam[0]])
        else:
            lt = vector([P_matrix[1,0]/lam[1], P_matrix[1,1]/lam[1]])
        return lam, lt
        
    # Otherwise, default to column method (or P|eta])
    c1 = vector(P_matrix.column(0))
    c2 = vector(P_matrix.column(1))
    
    if c1.norm() > 1e-10:
        lam = c1
        if abs(lam[0]) > 1e-10:
            lt = vector([P_matrix[0,0]/lam[0], P_matrix[0,1]/lam[0]])
        else:
            lt = vector([P_matrix[1,0]/lam[1], P_matrix[1,1]/lam[1]])
    elif c2.norm() > 1e-10:
        lam = c2
        if abs(lam[0]) > 1e-10:
            lt = vector([P_matrix[0,0]/lam[0], P_matrix[0,1]/lam[0]])
        else:
            lt = vector([P_matrix[1,0]/lam[1], P_matrix[1,1]/lam[1]])
    else:
        return vector([0,0]), vector([0,0])
        
    return lam, lt

if __name__ == "__main__":
    print("Testing BCFW implementation...")
    lambdas, tildes = sample_spinors_from_twistor(seed=42, n=6)
    
    # Test shift 0, 1
    # Channel s_012 is invariant under shift 0,1 because p0+p1 is invariant.
    # We need a channel where one index is in the set and one is out.
    # Try s_023 (0 is in, 1 is out).
    channel = [0, 2, 3]
    print(f"Testing shift on 0,1 for channel {channel}")
    
    z_star = solve_bcfw_pole(lambdas, tildes, 0, 1, channel)
    
    print(f"Found z_star: {z_star}")
    
    if z_star is not None:
        L_star, Lt_star = bcfw_shift_spinors(lambdas, tildes, 0, 1, z_star)
        s_val = get_channel_s(L_star, Lt_star, channel)
        print(f"s_channel(z_star) = {s_val}")
        
        if abs(s_val) < 1e-9:
            print("PASS: s_channel vanishes.")
        else:
            print("FAIL: s_channel did not vanish.")
    else:
        print("FAIL: z_star was None (constant channel?)")

