import sys
import os
from sage.all import *
import numpy as np

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.physics_limits.soft import compute_soft_ratio_n6

def extrapolate_soft_limit():
    print("M3.2: Extrapolating Soft Limit Convergence...")
    
    # epsilons = 10^-k
    ks = range(3, 9) # 3 to 8
    
    results = []
    
    print(f"{'epsilon':<15} | {'ratio':<20} | {'|1-ratio|':<20}")
    print("-" * 60)
    
    for k in ks:
        eps = QQ(1) / (10**k)
        ratio, status = compute_soft_ratio_n6(eps, seed=123)
        
        if status != "ok":
            print(f"{float(eps):<15.2e} | Error: {status}")
            continue
            
        diff = abs(1 - ratio)
        print(f"{float(eps):<15.2e} | {float(ratio):<20.10f} | {float(diff):<20.10e}")
        
        # Avoid log(0)
        if diff > 0:
            results.append((float(eps), float(diff)))
        
    # Fit diff ~ A * eps^p
    # log(diff) ~ log(A) + p * log(eps)
    
    if len(results) < 3:
        print("Not enough points to fit (or exact match).")
        # If all diffs are 0, that's perfect convergence
        return

    # Use least squares on log-log
    log_eps = np.log([r[0] for r in results])
    log_diff = np.log([r[1] for r in results])
    
    # y = mx + c => p * log_eps + log_A
    A_mat = np.vstack([log_eps, np.ones(len(log_eps))]).T
    p, c = np.linalg.lstsq(A_mat, log_diff, rcond=None)[0]
    
    print("-" * 60)
    print(f"Fit Result: |R - 1| ~ {np.exp(c):.4e} * epsilon^{p:.4f}")
    
    if p >= 0.9:
        print("PASS: Linear convergence (or better) detected.")
    else:
        print("FAIL: Convergence too slow (p < 0.9).")

if __name__ == "__main__":
    extrapolate_soft_limit()



