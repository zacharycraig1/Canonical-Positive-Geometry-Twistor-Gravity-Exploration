#!/usr/bin/env sage
# =============================================================================
# SPINOR HELICITY MODULE: Physical extraction from Momentum Twistors
# =============================================================================

from sage.all import *

def extract_spinors_from_twistor(twistor):
    """
    Extract spinor-helicity data from momentum twistor Z.
    Returns (lambdas, tildelambdas, x_points).
    """
    n = twistor.n
    Z = twistor.Z
    
    # 1. Extract Lambdas and Mus
    lambdas = []
    mus = []
    for i in range(n):
        lambdas.append(vector(QQ, [Z[i][0], Z[i][1]]))
        mus.append(vector(QQ, [Z[i][2], Z[i][3]]))
        
    # 2. Compute Region Momenta x_i
    # x_i satisfies mu_i = x_i lambda_i and mu_{i-1} = x_i lambda_{i-1}
    # Solution: x_i = (lambda_i mu_{i-1} - lambda_{i-1} mu_i) / <i-1 i>
    # Note: x_i is a 2x2 matrix (vector of 4 components? No, spinor indices)
    # x_i^{alpha beta} or x_i^{alpha dot_beta}?
    # In twistor space, Z = (lambda, mu). mu^dot_alpha = x^alpha dot_alpha lambda_alpha
    # Usually mu is dual spinor.
    # Let's use the formula:
    # x_i = ( |i> [mu_{i-1}| - |i-1> [mu_i| ) / <i-1 i>
    # where [mu| is just the mu vector treated as row?
    # Let's assume mu has index dot_alpha.
    
    x_points = []
    for i in range(n):
        im1 = (i - 1) % n
        denom = twistor.get_angle(im1, i)
        if denom == 0:
            return None, None, None
            
        lam_i = lambdas[i]
        lam_im1 = lambdas[im1]
        mu_i = mus[i]
        mu_im1 = mus[im1]
        
        # Outer products: lam * mu^T
        # x_i = (lam_i * mu_im1^T - lam_im1 * mu_i^T) / denom
        # Result is 2x2 matrix
        
        # Manual outer product
        m1 = matrix(QQ, 2, 2)
        m2 = matrix(QQ, 2, 2)
        for r in range(2):
            for c in range(2):
                m1[r,c] = lam_i[r] * mu_im1[c]
                m2[r,c] = lam_im1[r] * mu_i[c]
                
        x_i = (m1 - m2) / denom
        x_points.append(x_i)
        
    # 3. Compute Momenta p_i and Tilde Lambdas
    tildelambdas = []
    for i in range(n):
        ip1 = (i + 1) % n
        p_i_matrix = x_points[i] - x_points[ip1]
        
        # p_i = lambda_i * tilde_lambda_i^T
        # We know lambda_i. We need tilde_lambda_i.
        # p_i[r, c] = lam_i[r] * til_i[c]
        
        lam = lambdas[i]
        # Pick a non-zero component of lambda to solve
        if lam[0] != 0:
            til_0 = p_i_matrix[0, 0] / lam[0]
            til_1 = p_i_matrix[0, 1] / lam[0]
        elif lam[1] != 0:
            til_0 = p_i_matrix[1, 0] / lam[1]
            til_1 = p_i_matrix[1, 1] / lam[1]
        else:
            return None, None, None
            
        tildelambdas.append(vector(QQ, [til_0, til_1]))
        
    return lambdas, tildelambdas, x_points

def compute_physical_square_bracket(twistor, i, j):
    # This is now redundant if we have tildelambdas, but kept for compatibility
    # if called directly. Better to use extracted spinors.
    # But for now, let's just return None to force usage of extracted spinors
    return None

def compute_physical_mandelstam(twistor, i, j):
    # Redundant
    return None
