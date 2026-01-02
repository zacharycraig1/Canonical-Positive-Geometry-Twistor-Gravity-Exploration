import numpy as np
from scipy.optimize import root
import sys

def moment_map_and_jacobian(log_x, poly_coeffs, poly_exponents):
    """
    Computes X = grad_logx log p(x) and its Jacobian J = dX / dlogx.
    
    p(x) = sum c_v exp(v . log_x)
    """
    # log_x is vector of size d
    # poly_exponents is N x d matrix
    # poly_coeffs is size N
    
    # 1. Compute monomials and p(x)
    # log_m = V . log_x
    log_monomials = np.dot(poly_exponents, log_x)
    
    # Shift for stability? 
    # X doesn't change if we scale p(x) -> p(x)/M
    # So we can subtract max(log_monomials)
    shift = np.max(log_monomials)
    monomials_shifted = np.exp(log_monomials - shift)
    
    # p_shifted = sum c_v m_v
    p_val_shifted = np.dot(poly_coeffs, monomials_shifted)
    
    # weights w_v = c_v m_v / p_val
    # X = sum w_v v
    weights = (poly_coeffs * monomials_shifted) / p_val_shifted
    
    X = np.dot(weights, poly_exponents)
    
    # Jacobian J_ij = dX_i / d log x_j
    # X_i = sum_v w_v v_i
    # d w_v / d log x_j = w_v (v_j - X_j)
    # J_ij = sum_v v_i dw_v/dlogx_j = sum_v v_i w_v (v_j - X_j)
    #      = sum_v w_v v_i v_j - (sum w_v v_i) X_j
    #      = <v_i v_j> - <v_i> <v_j>
    #      = Covariance of v under measure w
    
    # Covariance matrix
    # E[v v^T] - E[v] E[v]^T
    # E[v] = X
    
    # Weighted sum of outer products
    # V is N x d. weights is N.
    # We want V^T W V  where W is diag(weights)
    
    # Efficient computation:
    # J = (V.T * weights) @ V - outer(X, X)
    
    weighted_V = poly_exponents.T * weights # broadcast weights
    second_moment = np.dot(weighted_V, poly_exponents)
    
    J = second_moment - np.outer(X, X)
    
    return X, J

def solve_saddle(X_target, poly_coeffs, poly_exponents, initial_guess=None):
    """
    Solves X(log_x) = X_target for log_x.
    Returns solution log_x_star.
    """
    dim = len(X_target)
    if initial_guess is None:
        initial_guess = np.zeros(dim)
        
    def func(log_x):
        X, _ = moment_map_and_jacobian(log_x, poly_coeffs, poly_exponents)
        return X - X_target
        
    def jac(log_x):
        _, J = moment_map_and_jacobian(log_x, poly_coeffs, poly_exponents)
        return J
        
    # Try different methods if hybr fails
    methods = ['hybr', 'lm', 'broyden1']
    
    for method in methods:
        try:
            if method == 'lm':
                # lm requires least_squares or specialized call, root(method='lm') works for N=N
                # But jacobian signature might differ or be required
                sol = root(func, initial_guess, jac=jac, method='lm', tol=1e-9)
            else:
                sol = root(func, initial_guess, jac=jac, method=method, tol=1e-9)
            
            if sol.success:
                return sol.x, True
        except Exception:
            continue
            
    # If standard attempts fail, try random restarts with hybr
    for i in range(10):
        # Scale guess?
        guess = np.random.randn(dim) * 2.0
        try:
            sol = root(func, guess, jac=jac, method='hybr', tol=1e-9)
            if sol.success:
                return sol.x, True
        except Exception:
            continue

    return initial_guess, False

def compute_pushforward_saddle(X_target, poly_coeffs, poly_exponents):
    """
    Computes sum 1 / det(J(x*)) over solutions.
    For moment map of positive polynomial, solution is unique.
    
    Returns:
        float: Value of 1/det(J) at saddle point.
    """
    poly_exponents = np.array(poly_exponents)
    poly_coeffs = np.array(poly_coeffs)
    X_target = np.array(X_target)
    
    log_x_star, success = solve_saddle(X_target, poly_coeffs, poly_exponents)
    
    if not success:
        # Try a few random restarts if failed?
        for i in range(5):
             guess = np.random.randn(len(X_target))
             log_x_star, success = solve_saddle(X_target, poly_coeffs, poly_exponents, initial_guess=guess)
             if success: break
             
    if not success:
        raise ValueError("Could not solve saddle point equation.")
        
    _, J = moment_map_and_jacobian(log_x_star, poly_coeffs, poly_exponents)
    
    det_J = np.linalg.det(J)
    
    # Canonical form is positive? 
    # J is positive definite for moment map (covariance matrix).
    # det(J) > 0.
    
    return 1.0 / det_J


