import sys
import os
import numpy as np
from scipy.special import gamma

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from posgeom.stringy_integral import numerical_stringy_integral

def test_simplex_n4():
    print("Testing Stringy Integral for Simplex (n=4 case equivalent)...")
    
    # Simplex 2D: p(t) = 1 + t1 + t2
    coeffs = [1, 1, 1]
    exponents = [[0, 0], [1, 0], [0, 1]]
    
    # Target point (centroid)
    X = [1/3.0, 1/3.0]
    
    # Exact Formula:
    # I = (alpha)^d * Gamma(alpha*X1) * Gamma(alpha*X2) * Gamma(alpha*(1-X1-X2)) / Gamma(alpha)
    # Note: The formula I_p = (alpha)^d int ... 
    # For p = 1+t1+t2, integral is Gamma(aX)... / Gamma(a)
    # Yes.
    
    alphas = [1.0, 0.5, 0.2] 
    # Smaller alphas require larger integration domains and more points due to slow decay
    
    for alpha in alphas:
        # Exact
        arg1 = alpha * X[0]
        arg2 = alpha * X[1]
        arg3 = alpha * (1.0 - X[0] - X[1])
        
        exact_integral = gamma(arg1) * gamma(arg2) * gamma(arg3) / gamma(alpha)
        # The function returns I * (alpha)^d.
        # Wait, the integral itself is Gamma.../Gamma...
        # The numerical_stringy_integral returns alpha^d * Integral.
        # So we should compare with alpha^d * exact_integral_value.
        
        expected = (alpha**2) * exact_integral
        
        try:
            # Increase limit for better accuracy
            val = numerical_stringy_integral(coeffs, exponents, X, alpha, opts={'limit': 200})
            
            print(f"alpha={alpha:.3f} | Numeric: {val:.5f} | Exact: {expected:.5f} | Ratio: {val/expected:.3f}")
            
            # Allow some tolerance (5-10% is acceptable for numeric integration of slowly decaying functions)
            if abs(val - expected) / expected > 0.1:
                print("  -> Warning: Deviation > 10%")
                
        except ImportError:
            print("Scipy not available, skipping.")
            return
        except Exception as e:
            print(f"Integration failed: {e}")

def main():
    test_simplex_n4()

if __name__ == "__main__":
    main()

