import sys
import os
import json
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_polynomial
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian

def phaseG4_check_polynomial_value():
    print("Phase G4: Verifying Forest Polynomial on Physical Kinematics...")
    
    n = 6
    roots = [0, 1, 2]
    
    # 1. Get Polynomial
    print("Generating Forest Polynomial (symbolic)...")
    F_poly = get_forest_polynomial(n, roots)
    print(f"  Polynomial has {len(F_poly.monomials())} terms.")
    
    # 2. Check on seeds
    n_seeds = 3
    for seed in range(n_seeds):
        print(f"\nProcessing seed {seed}...")
        try:
            lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=n)
        except: continue
        
        # References
        x = vector(QQ, [1, 0])
        y = vector(QQ, [0, 1])
        
        # A) Evaluate z_ij
        z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        
        # B) Evaluate Polynomial F(z)
        # Construct substitution dict using ring generators
        R = F_poly.parent()
        gens = R.gens_dict()
        subs = {}
        for (i, j), val in z_map.items():
            key = f"z_{i}_{j}"
            if key in gens:
                subs[gens[key]] = val
            
        F_val = F_poly.subs(subs)
        
        # If result is still a polynomial (constant), cast it
        if hasattr(F_val, "constant_coefficient"):
            F_val = F_val.constant_coefficient()
            
        # C) Compare with Determinant
        L_tilde, _, _ = hodges_weighted_laplacian(lambdas, tildes, x, y)
        indices = [3, 4, 5]
        det_val = L_tilde.matrix_from_rows_and_columns(indices, indices).det()
        
        # Check equality
        # Note: L_ij = -z_ij.
        # In the determinant expansion of a 3x3 matrix, cubic terms involve 3 entries.
        # Product of (-z) is (-1)^3 z^3 = -z^3.
        # But wait, expansion of det(L) also has signs from permutation.
        # Matrix Tree Theorem says: det(L^(R)) = Sum_{Forests} Prod_{edges} w_e.
        # Here w_e in L is z_e (with minus sign).
        # Standard result: If L = D - A, then det = Sum weights.
        # Our L has -z_ij on off-diagonals. So A_ij = z_ij.
        # So it should be EXACTLY F_val.
        
        diff = abs(F_val - det_val)
        print(f"  Poly Val: {F_val}")
        print(f"  Det Val:  {det_val}")
        
        if diff < 1e-9:
            print("  MATCH")
        else:
            print(f"  MISMATCH (Diff: {diff})")

if __name__ == "__main__":
    phaseG4_check_polynomial_value()

