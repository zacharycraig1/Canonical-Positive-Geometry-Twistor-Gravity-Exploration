import numpy as np
import sys

# Try to import scipy
try:
    from scipy.integrate import nquad
    has_scipy = True
except ImportError:
    has_scipy = False
    sys.stderr.write("Warning: scipy not found, stringy integral will fail.\n")

def numerical_stringy_integral(poly_coeffs, poly_exponents, X_target, alpha_prime, opts=None):
    """
    Computes the stringy integral I = (alpha')^d * int_{R>0^d} prod(dx/x) x^(alpha' X) p(x)^(-alpha').
    
    Args:
        poly_coeffs (list): Coefficients c_v of polynomial p(x) = sum c_v x^v
        poly_exponents (list of list): Exponent vectors v. Must be dimension d.
        X_target (list): Target point X in intrinsic space (dim d).
        alpha_prime (float): String tension parameter.
        opts (dict): Options for integration.
        
    Returns:
        float: The integral value.
    """
    if not has_scipy:
        raise ImportError("scipy.integrate is required for numerical_stringy_integral")

    dim = len(X_target)
    poly_exponents = np.array(poly_exponents)
    poly_coeffs = np.array(poly_coeffs)
    X_target = np.array(X_target)
    
    # Pre-check dimensions
    if poly_exponents.shape[1] != dim:
        raise ValueError(f"Exponents dimension {poly_exponents.shape[1]} != target dimension {dim}")

    def integrand(*args):
        # args is tuple (t_1, ..., t_d)
        t = np.array(args)
        
        # Avoid zero or negative (though nquad shouldn't pass them if bounds are (0, inf))
        if np.any(t <= 1e-12): 
            # Regularize slightly near 0 if needed, or return 0 if measures zero
            # But t^(negative) blows up.
            # nquad usually handles open interval (0, inf) by transformation?
            # We assume t > 0.
            return 0.0
            
        # Term 1: t^(alpha' X - 1)
        # We compute prod(t_i ^ (alpha * X_i - 1))
        # = exp( sum (alpha X_i - 1) log t_i )
        log_t = np.log(t)
        term1_log = np.sum((alpha_prime * X_target - 1.0) * log_t)
        
        # Polynomial value p(t)
        # sum c_v exp(v . log t)
        # Compute v . log t for all v
        log_monomials = np.dot(poly_exponents, log_t)
        # Shift to avoid overflow/underflow in exp if needed?
        # p(t) is sum of positive terms usually (coeffs > 0 for pos geom).
        # Assuming positive coefficients for now (or at least p(t) > 0).
        monomials = np.exp(log_monomials)
        p_val = np.dot(poly_coeffs, monomials)
        
        if p_val <= 0:
            # Should not happen for positive geometry polynomials in positive octant
            return 0.0
            
        # Term 2: p(t)^(-alpha') = exp(-alpha' * log p_val)
        term2_log = -alpha_prime * np.log(p_val)
        
        total_log = term1_log + term2_log
        
        # Apply prefactor later to avoid scaling issues?
        # Function returns integrand value.
        return np.exp(total_log)

    # Scale result by (alpha')^d
    # Bounds (0, inf)
    ranges = [(0, np.inf)] * dim
    
    # Options
    # Increasing limit helps for slowly decaying functions
    run_opts = {'limit': 100}
    if opts:
        run_opts.update(opts)
        
    val, error = nquad(integrand, ranges, opts=run_opts)
    
    return val * (alpha_prime ** dim)

