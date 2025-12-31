import sys
import os
import json
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian
from src.chy_oracle.forest_enumerator import enumerate_rooted_forests, get_forest_polynomial
from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket

def verify_forest_expansion():
    print("Starting Phase F2.1: Forest Expansion Verification...")
    
    n_seeds = 3
    results_log = []
    
    # Precompute forests for n=6, roots={0,1,2}
    # This is purely combinatorial
    print("Enumerating forests...")
    roots = [0, 1, 2]
    forests = enumerate_rooted_forests(6, roots)
    print(f"Found {len(forests)} rooted forests.")
    
    for seed in range(n_seeds):
        print(f"\nProcessing seed {seed}...")
        
        try:
            lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=6)
        except Exception:
            continue
            
        # Reference spinors
        x = vector(QQ, [1, 0])
        y = vector(QQ, [0, 1])
        
        try:
            L_tilde, C, _ = hodges_weighted_laplacian(lambdas, tildes, x, y)
        except ValueError:
            continue
            
        # Recompute weights for validation
        weights_dict = {}
        for i in range(6):
            for j in range(i+1, 6):
                ang = ang_bracket(lambdas[i], lambdas[j])
                sq = sq_bracket(tildes[i], tildes[j])
                if ang != 0:
                    weights_dict[(i, j)] = sq / ang

        # 1. Compute Determinant Minor
        indices = [3, 4, 5]
        det_val = L_tilde.matrix_from_rows_and_columns(indices, indices).det()
        
        # 2. Compute Forest Sum
        # weights_dict is (i,j) -> w_ij
        # C is list
        
        def w_func(u, v):
            # Sort u, v because weights_dict key might be sorted
            if u > v: u, v = v, u
            return weights_dict.get((u, v), 0)
            
        def c_func(u):
            return C[u]
            
        forest_sum = get_forest_polynomial(forests, w_func, c_func)
        
        # Compare
        # The determinant expansion of L_{sub} where L_ij = -wCC is sum(prod(wCC)).
        # Because (-1)^(n-k) matches? 
        # Size of minor is 3x3. Expansion has terms of degree 3.
        # Prod(-wCC) = (-1)^3 * prod(wCC) = - prod(wCC).
        # But wait, diagonal L_ii = sum w_ik C_i C_k (positive terms).
        # Determinant of Laplacian minor is POSITIVE sum of forests?
        # Standard Matrix Tree Theorem: det(L^(R)) = sum of forests.
        # With signs: Off-diagonals are negative. 
        # If L is standard Laplacian (M-matrix), det is positive.
        # Our L has off-diagonals -wCC. w is scalar (could be negative but we found positive on moment curve).
        # C^2 is positive? No, C is complex. C^2 is complex.
        # But we compare algebraic equality, so signs must match exactly.
        
        diff = det_val - forest_sum
        ratio = det_val / forest_sum if forest_sum != 0 else 0
        
        print(f"  Det: {det_val}")
        print(f"  Sum: {forest_sum}")
        
        if abs(diff) < 1e-9:
            print("  MATCH")
            status = "MATCH"
        elif abs(det_val + forest_sum) < 1e-9:
             print("  MATCH (Sign Flip)")
             status = "MATCH_SIGN_FLIP"
        else:
             print(f"  MISMATCH (Ratio {float(abs(ratio))})")
             status = "MISMATCH"
             
        results_log.append({
            "seed": seed,
            "status": status,
            "det": str(det_val),
            "sum": str(forest_sum)
        })

    with open("results/phaseF_forest_vs_det.json", "w") as f:
        json.dump(results_log, f, indent=2)

if __name__ == "__main__":
    verify_forest_expansion()

