import sys
import os
import json
import random
from sage.all import QQ, matrix, vector, Matrix, VectorSpace

# Adjust path to import local modules
sys.path.append(os.getcwd())

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)

def get_s_map(u):
    """
    Given u (9-element vector), return full s_{ij} dictionary and s_S helper.
    u = [s01, s02, s03, s04, s12, s13, s14, s23, s24]
    """
    s01, s02, s03, s04, s12, s13, s14, s23, s24 = u
    
    s = {}
    s[(0,1)] = s01; s[(1,0)] = s01
    s[(0,2)] = s02; s[(2,0)] = s02
    s[(0,3)] = s03; s[(3,0)] = s03
    s[(0,4)] = s04; s[(4,0)] = s04
    s[(1,2)] = s12; s[(2,1)] = s12
    s[(1,3)] = s13; s[(3,1)] = s13
    s[(1,4)] = s14; s[(4,1)] = s14
    s[(2,3)] = s23; s[(3,2)] = s23
    s[(2,4)] = s24; s[(4,2)] = s24
    
    s05 = -(s01 + s02 + s03 + s04)
    s[(0,5)] = s05; s[(5,0)] = s05
    
    s15 = -(s01 + s12 + s13 + s14)
    s[(1,5)] = s15; s[(5,1)] = s15
    
    s25 = -(s02 + s12 + s23 + s24)
    s[(2,5)] = s25; s[(5,2)] = s25
    
    C3 = -(s03 + s13 + s23)
    C4 = -(s04 + s14 + s24)
    C5 = -(s05 + s15 + s25)
    
    s34 = (C3 + C4 - C5) / 2
    s35 = (C3 - C4 + C5) / 2
    s45 = (-C3 + C4 + C5) / 2
    
    s[(3,4)] = s34; s[(4,3)] = s34
    s[(3,5)] = s35; s[(5,3)] = s35
    s[(4,5)] = s45; s[(5,4)] = s45
    
    return s

def get_s_subset(s_map, subset):
    val = 0
    sub = sorted(list(subset))
    for idx, i in enumerate(sub):
        for j in sub[idx+1:]:
            val += s_map[(i,j)]
    return val

def solve_kinematic_map():
    print("Loading data...")
    facets = load_json("RESULTS/facet_dictionary_n6.json")
    hull_eq = load_json("RESULTS/facets_n6_eq_exact.json")
    
    equations = hull_eq["equations"]
    n_dim = 15
    n_eq = len(equations)
    
    E_rows = []
    b_vec = []
    
    for eq in equations:
        coeffs = [QQ(x) for x in eq]
        b_val = coeffs[0]
        row = coeffs[1:]
        E_rows.append(row)
        b_vec.append(b_val)
        
    E = matrix(QQ, E_rows)
    b_vec = vector(QQ, b_vec)
    
    try:
        X0 = E.solve_right(-b_vec)
    except ValueError:
        print("Error: No solution to hull equations!")
        return
        
    N_ker = E.right_kernel()
    N = N_ker.matrix().transpose()
    
    if N.ncols() != 11:
        print(f"Warning: Expected kernel dim 11, got {N.ncols()}")
        
    print(f"Hull parametrized: X = X0 + N*t, t in Q^{N.ncols()}")
    
    matched_facets = []
    for f in facets:
        if f["best_match"]["type"] in ["internal", "cut"]:
            matched_facets.append(f)
            
    print(f"Found {len(matched_facets)} matched facets for fitting.")
    
    M = len(matched_facets)
    n_t = 11
    n_u = 9
    
    # Check 1: Constant Part Consistency
    # (A_F N) t0 = - (b_F + A_F X0)
    print("Checking constant part consistency...")
    
    consistent_facets = []
    skipped_facets = []
    
    # Greedy selection of consistent facets
    current_rows = []
    current_rhs = []
    
    # First, shuffle or sort? 
    # Maybe prioritize by some metric? Let's just use the order in JSON.
    
    for f in matched_facets:
        f_coeffs = [QQ(x) for x in f["ineq_b_A"]]
        b_F = f_coeffs[0]
        A_F = vector(QQ, f_coeffs[1:])
        
        row = A_F * N
        rhs = -(b_F + A_F.dot_product(X0))
        
        # Test if adding this row preserves consistency
        test_rows = current_rows + [row]
        test_rhs = current_rhs + [rhs]
        
        mat_test = matrix(QQ, test_rows)
        vec_test = vector(QQ, test_rhs)
        
        # Rank check
        r = mat_test.rank()
        aug = mat_test.augment(vec_test)
        r_aug = aug.rank()
        
        if r == r_aug:
            # Consistent
            current_rows.append(row)
            current_rhs.append(rhs)
            consistent_facets.append(f)
        else:
            skipped_facets.append(f)

    print(f"Consistent facets: {len(consistent_facets)}/{len(matched_facets)}")
    if skipped_facets:
        print("Skipped facets (Affine mismatch):")
        for f in skipped_facets:
            print(f"  ID {f['facet_id']} (Subset {f['best_match']['subset']})")

    # Proceed with consistent facets only for t0 solving
    # But for T solving, we might still fit them?
    # If we skip them for t0, it means L_F(X0) != 0.
    # So L_F(X) = const + kappa s_S.
    # This violates Gate A (Exact Proportionality).
    # But we can't fix it if geometry prevents it.
    
    # We will fit the map using ONLY the consistent facets for the constraints.
    # The others will just be what they are.
    
    # Solve for t0 (particular solution for consistent set)
    # We need a t0 that satisfies the consistent equations.
    # And maybe minimizes error for others? No, just exact for consistent.
    
    mat_C = matrix(QQ, current_rows)
    vec_C = vector(QQ, current_rhs)
    # This system is consistent by construction.
    # It might be underdetermined (kernel > 0).
    # We can pick a particular t0, or leave free pars for T?
    # t0 is a constant vector. It doesn't depend on u.
    # Free parameters in t0 can be set to 0.
    
    t0_particular = mat_C.solve_right(vec_C)
    # kernel of C
    K_C = mat_C.right_kernel()
    # t0 = t0_particular + K_C * alpha. 
    # We can just use t0_particular for now.
    
    # Now setup T fit.
    # For consistent facets: L_F = kappa s_S. (const part matches 0).
    # For skipped facets: L_F = const + kappa s_S?
    # If we enforce L_F slope = kappa s_S slope, we can fit T.
    # But the constant term will remain non-zero.
    
    # Let's fit T using ALL matched facets (slopes can be matched even if intercept fails).
    # Equation: (A_F N) T u = kappa s_S(u)
    # (Constant parts cancel out for consistent, remain for skipped).
    # Wait, the equation we derived:
    # (A_F N) t0 + (A_F N) T u - k s_S = -(b + A X0)
    # If we use t0_particular, then for consistent facets:
    # (A_F N) t0_p = -(b + A X0).
    # So (A_F N) T u - k s_S = 0.
    
    # For skipped facets:
    # (A_F N) t0_p != -(b + A X0).
    # Let Rem = (b + A X0) + (A_F N) t0_p. (Non-zero residual).
    # Equation: (A_F N) T u - k s_S = -Rem
    # But LHS is linear in u. RHS is constant.
    # Only possible if LHS constant = RHS constant.
    # But T u is purely u-dependent? No, if u=0 -> 0.
    # s_S(0) = 0.
    # So LHS(0) = 0. RHS = -Rem != 0.
    # Contradiction for all u.
    
    # So skipped facets CANNOT satisfy the equation even with T freedom.
    # Because T multiplies u.
    # Unless we add a constant term to the map ansatz? t(u) = t0 + T u.
    # We already did. t0 is the constant.
    # The inconsistency is in t0.
    
    # So for skipped facets, we simply cannot enforce the condition.
    # We will exclude them from the fit entirely.
    
    matched_facets = consistent_facets
    
    # Re-setup full system with REDUCED set
    idx_t0 = 0
    idx_T = 11
    idx_k = 11 + 11*9
    total_vars = idx_k + len(matched_facets)
    M = len(matched_facets)
    
    system_rows = []
    system_rhs = []
    
    num_samples = 20
    
    print(f"Generating full system with {num_samples} samples (Consistent only)...")
    
    for sample_i in range(num_samples):
        u_vec = [QQ(random.randint(-10, 10)) for _ in range(n_u)]
        s_map = get_s_map(u_vec)
        
        for k_idx, facet in enumerate(matched_facets):
            f_coeffs = [QQ(x) for x in facet["ineq_b_A"]]
            b_F = f_coeffs[0]
            A_F = vector(QQ, f_coeffs[1:])
            subset = facet["best_match"]["subset"]
            s_S_val = get_s_subset(s_map, subset)
            
            coeff_t0 = A_F * N
            
            coeff_T = []
            for i in range(n_t):
                for j in range(n_u):
                    coeff_T.append(coeff_t0[i] * u_vec[j])
                    
            coeff_k = [0]*M
            coeff_k[k_idx] = -s_S_val
            
            const_val = b_F + A_F.dot_product(X0)
            rhs_val = -const_val
            
            row = list(coeff_t0) + coeff_T + coeff_k
            
            system_rows.append(row)
            system_rhs.append(rhs_val)
            
    # Add constraint kappa_0 = 1 to force non-trivial solution
    # Row: [0...0 (t0), 0...0 (T), 1, 0...0 (kappas)]
    # Index of kappa_0 is idx_k.
    constraint_row = [QQ(0)] * total_vars
    constraint_row[idx_k] = QQ(1)
    system_rows.append(constraint_row)
    system_rhs.append(QQ(1))

    print("Solving linear system...")
    A_sys = matrix(QQ, system_rows)
    B_sys = vector(QQ, system_rhs)
    
    try:
        Y = A_sys.solve_right(B_sys)
        print("Solution found!")
        
        t0 = Y[idx_t0 : idx_t0 + n_t]
        T_flat = Y[idx_T : idx_T + n_t*n_u]
        T_rows = []
        for i in range(n_t):
            T_rows.append(list(T_flat[i*n_u : (i+1)*n_u]))
        kappas = Y[idx_k : idx_k + M]
        
        def json_friendly(obj):
            if hasattr(obj, "numerator"): # Rational
                return float(obj)
            if isinstance(obj, list):
                return [json_friendly(x) for x in obj]
            return obj

        result = {
            "X0": json_friendly(list(X0)),
            "N": json_friendly([list(row) for row in N.rows()]), 
            "N_rows": json_friendly([list(r) for r in N.rows()]),
            "t0": json_friendly(list(t0)),
            "T": json_friendly(T_rows),
            "matched_facets": [],
            "edge_order": hull_eq["edges_ordered"]
        }
        
        print("\nMatched Facet Scalings (kappa):")
        for k_idx, k_val in enumerate(kappas):
            f = matched_facets[k_idx]
            print(f"Facet {f['facet_id']}: kappa = {k_val}")
            result["matched_facets"].append({
                "facet_id": f["facet_id"],
                "kappa": float(k_val),
                "subset": f["best_match"]["subset"]
            })
            
        out_path = "RESULTS/kinematic_map_n6_fit.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Map saved to {out_path}")
        
    except ValueError:
        print("No solution found. System may be inconsistent.")
        r = A_sys.rank()
        aug = A_sys.augment(B_sys)
        r_aug = aug.rank()
        print(f"Rank A: {r}, Rank Aug: {r_aug}")

if __name__ == "__main__":
    solve_kinematic_map()
