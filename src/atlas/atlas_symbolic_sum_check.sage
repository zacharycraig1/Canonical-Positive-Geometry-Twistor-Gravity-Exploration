import sys
import os
import json
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

def run_symbolic_check():
    print("Running Algebraic Sum Check (Numeric Evaluation of Polynomials)...")
    
    # Load Coeffs
    coeff_path = "RESULTS/atlas_coeffs_exact.json"
    if not os.path.exists(coeff_path):
        print("Coefficients not found.")
        return
        
    with open(coeff_path, 'r') as f:
        data = json.load(f)
        
    coeffs = {}
    for item in data:
        roots = tuple(sorted(item['roots']))
        coeffs[roots] = item['coeff']
        
    n = 6
    RF = RealField(200)
    
    # Define Evaluation Function
    def evaluate_det_L_R(roots, z_map):
        # Build L (n x n)
        # L_ij = -z_ij
        L = matrix(RF, n, n)
        for i in range(n):
            diag = 0
            for j in range(n):
                if i == j: continue
                uv = tuple(sorted((i,j)))
                val = z_map[uv]
                L[i,j] = -val
                diag += val
            L[i,i] = diag
            
        non_roots = sorted([i for i in range(n) if i not in roots])
        L_R = L[non_roots, non_roots]
        return L_R.det()

    # Reference Roots (Pick the first one with non-zero coeff)
    ref_roots = None
    for r, c in coeffs.items():
        if abs(c) > 0.1:
            ref_roots = r
            break
            
    if ref_roots is None:
        print("No non-zero coefficients found.")
        return
        
    print(f"Reference Roots: {ref_roots}")
    
    print("\nTesting on random positive z (Schwarz-Zippel check)...")
    
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i,j))
            
    num_trials = 10
    ratios = []
    
    for k in range(num_trials):
        z_map = {uv: RF(rnd.uniform(0.1, 2.0)) for uv in edges}
        
        # LHS: Sum c_R * det(L_R)
        lhs = RF(0)
        for roots, c in coeffs.items():
            if c == 0: continue
            det_val = evaluate_det_L_R(roots, z_map)
            lhs += c * det_val
            
        # RHS: det(L_ref)
        rhs = evaluate_det_L_R(ref_roots, z_map)
        
        ratio = lhs / rhs if abs(rhs) > 1e-20 else 0
        print(f"Trial {k}: LHS={lhs:.2e}, RHS={rhs:.2e}, Ratio={ratio:.4f}")
        ratios.append(ratio)
        
    avg = sum(ratios) / len(ratios)
    dev = sum([(r - avg)**2 for r in ratios]) / len(ratios)
    
    print(f"\nAverage Ratio: {avg:.6f}")
    print(f"Variance: {dev:.2e}")
    
    result = {
        "status": "PASS" if dev < 1e-10 else "FAIL",
        "ratio": float(avg),
        "variance": float(dev),
        "ref_roots": list(ref_roots)
    }
    
    with open("RESULTS/algebraic_check_result.json", "w") as f:
        json.dump(result, f, indent=2)

if __name__ == "__main__":
    run_symbolic_check()


