import json
import os
import sys
import itertools
import math

def gcd_list(lst):
    """Compute GCD of a list of integers."""
    if not lst:
        return 1
    result = abs(lst[0])
    for x in lst[1:]:
        result = math.gcd(result, abs(x))
    return result

def normalize_vector(vec):
    """
    Normalize integer vector by dividing by GCD.
    Returns (normalized_vec, factor).
    Handles zero vector.
    """
    if all(v == 0 for v in vec):
        return tuple(vec), 1
    
    g = gcd_list(vec)
    return tuple(v // g for v in vec), g

def vector_sub(v1, v2):
    return [a - b for a, b in zip(v1, v2)]

def generate_subset_vectors(n, edges_ordered):
    """
    Generates cut and internal vectors for all non-trivial subsets.
    Returns a dict: vector_tuple -> list of (type, subset)
    """
    vectors = {}
    
    edge_to_idx = {tuple(e): i for i, e in enumerate(edges_ordered)}
    num_edges = len(edges_ordered)
    
    for r in range(1, n): 
        for subset in itertools.combinations(range(n), r):
            S = set(subset)
            S_list = sorted(list(subset))
            
            # 1. Cut Vector
            cut_vec = [0] * num_edges
            for u, v in edges_ordered:
                idx = edge_to_idx[(u, v)]
                if (u in S) != (v in S):
                    cut_vec[idx] = 1
            
            norm_cut, _ = normalize_vector(cut_vec)
            if norm_cut not in vectors:
                vectors[norm_cut] = []
            vectors[norm_cut].append({"type": "cut", "subset": S_list})

            # 2. Internal Vector
            if len(S) > 1:
                int_vec = [0] * num_edges
                for u, v in edges_ordered:
                    idx = edge_to_idx[(u, v)]
                    if u in S and v in S:
                        int_vec[idx] = 1
                
                norm_int, _ = normalize_vector(int_vec)
                if norm_int not in vectors:
                    vectors[norm_int] = []
                vectors[norm_int].append({"type": "internal", "subset": S_list})

    return vectors

def identify_facets():
    input_path = "RESULTS/facets_n6_ineq_exact.json"
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    print(f"Loading facets from {input_path}...")
    with open(input_path, "r") as f:
        data = json.load(f)
        
    n = data["n"]
    edges_ordered = data["edges_ordered"]
    roots = set(data.get("roots", [0, 1, 2]))
    facets = data["inequalities"] 
    
    print(f"Generating reference vectors for n={n}...")
    ref_vectors = generate_subset_vectors(n, edges_ordered)
    
    # Basis vectors
    basis_vectors = {}
    for i, edge in enumerate(edges_ordered):
        vec = [0]*len(edges_ordered)
        vec[i] = 1
        basis_vectors[tuple(vec)] = edge

    # Identify root-root edges
    rr_indices = set()
    for idx, (u, v) in enumerate(edges_ordered):
        if u in roots and v in roots:
            rr_indices.add(idx)
            
    print(f"Ignoring root-root edge indices: {sorted(list(rr_indices))}")

    # Helper to mask vector
    def mask_rr(vec):
        v = list(vec)
        for idx in rr_indices:
            v[idx] = 0
        return tuple(v)

    # Re-process reference vectors with masking
    masked_ref_vectors = {}
    for vec, candidates in ref_vectors.items():
        masked_vec, _ = normalize_vector(mask_rr(vec))
        if masked_vec not in masked_ref_vectors:
            masked_ref_vectors[masked_vec] = []
        masked_ref_vectors[masked_vec].extend(candidates)
        
    # Re-process basis vectors with masking
    masked_basis_vectors = {}
    for vec, edge in basis_vectors.items():
        masked_vec, _ = normalize_vector(mask_rr(vec))
        if masked_vec not in masked_basis_vectors:
            masked_basis_vectors[masked_vec] = edge

    # Global Sum Vector (masked)
    global_sum_vec = [0] * len(edges_ordered)
    for i in range(len(edges_ordered)):
        if i not in rr_indices:
            global_sum_vec[i] = 1
    # Note: sum(x) = 3 for n=6
    global_rhs = 3

    matches = []
    
    for i, f_coeffs_str in enumerate(facets):
        f_coeffs = [int(x) for x in f_coeffs_str]
        b = f_coeffs[0]
        A = list(f_coeffs[1:])
        
        # Zero out root-root coefficients in A
        for idx in rr_indices:
            A[idx] = 0
            
        norm_A, g_A = normalize_vector(A)
        neg_norm_A, g_neg_A = normalize_vector([-x for x in A])
        
        match_data = {
            "facet_id": i,
            "ineq": f_coeffs, 
            "identification": "unknown",
            "proof": None
        }
        
        # Case 1: Lower Bound x_e >= 0
        if norm_A in masked_basis_vectors:
            edge = masked_basis_vectors[norm_A]
            if b == 0:
                match_data["identification"] = "lower_bound"
                match_data["proof"] = {"edge": edge, "bound": 0}
        
        # Case 2: Upper Bound x_e <= 1
        elif neg_norm_A in masked_basis_vectors:
            edge = masked_basis_vectors[neg_norm_A]
            if b == 1:
                match_data["identification"] = "upper_bound"
                match_data["proof"] = {"edge": edge, "bound": 1}

        # Case 3: Subset Upper Bound (Internal or Cut)
        elif neg_norm_A in masked_ref_vectors:
            candidates = masked_ref_vectors[neg_norm_A]
            best_cand = candidates[0]
            # Try to refine best_cand
            for cand in candidates:
                S = set(cand["subset"])
                r_cap_s = len(S.intersection(roots))
                predicted_rhs = len(S) - max(1, r_cap_s)
                if cand["type"] == "internal" and b == predicted_rhs:
                     best_cand = cand
                     best_cand["predicted_rhs"] = predicted_rhs
                     break
            
            match_data["identification"] = f"subset_{best_cand['type']}_upper"
            match_data["proof"] = {
                "subset": best_cand["subset"],
                "rhs": b,
                "all_candidates": candidates
            }

        # Case 4: Subset Lower Bound
        elif norm_A in masked_ref_vectors:
            candidates = masked_ref_vectors[norm_A]
            best = candidates[0] 
            match_data["identification"] = f"subset_{best['type']}_lower"
            match_data["proof"] = {
                "subset": best["subset"],
                "rhs": -b,
                "all_candidates": candidates
            }

        # Case 5: Complement Logic (Implicit via Global Sum)
        if match_data["identification"] == "unknown":
            # Check Complement of A: C = Ones - A (if A is 0/1)
            # Or C = Ones + A (if A is 0/-1)
            
            # If A is positive (lower bound like), check if C = Ones - A is in Ref/Basis
            # sum(A) >= -b <=> sum(All) - sum(C) >= -b <=> 3 - sum(C) >= -b <=> sum(C) <= 3 + b
            # This converts Lower Bound on A to Upper Bound on C.
            
            # Construct candidate Complement vectors
            # Try C = Ones - A
            C_vec = vector_sub(global_sum_vec, A)
            norm_C, _ = normalize_vector(mask_rr(C_vec))
            
            if norm_C in masked_basis_vectors:
                 # sum(C) <= 3 + b
                 # If C is basis e_k, then x_k <= 3 + b
                 edge = masked_basis_vectors[norm_C]
                 match_data["identification"] = "implicit_upper_bound"
                 match_data["proof"] = {"edge": edge, "bound": 3+b}

            elif norm_C in masked_ref_vectors:
                 # sum(C) <= 3 + b
                 candidates = masked_ref_vectors[norm_C]
                 match_data["identification"] = "implicit_subset_upper"
                 match_data["proof"] = {"subset": candidates[0]["subset"], "rhs": 3+b}
                 
            # Try C = Ones + A (if A is negative, upper bound like)
            # sum(A) >= -b. A is negative. Let A = -A'. sum(-A') >= -b => sum(A') <= b.
            # If A' is Complement of C.
            # Effectively, A corresponds to "All - C" with -1 coeffs?
            # A = - (Ones - C) = C - Ones.
            # sum(C - Ones) >= -b => sum(C) - 3 >= -b => sum(C) >= 3 - b.
            # This converts Upper Bound on A (technically A is negative so sum(abs(A)) <= b)
            # to Lower Bound on C.
            
            # Check if A = C - Ones
            C_vec2 = vector_sub(A, [-x for x in global_sum_vec]) # A + Ones
            norm_C2, _ = normalize_vector(mask_rr(C_vec2))
            
            if norm_C2 in masked_basis_vectors:
                 # sum(C) >= 3 - b
                 edge = masked_basis_vectors[norm_C2]
                 # If 3-b = 0, it's x_e >= 0
                 match_data["identification"] = "implicit_lower_bound"
                 match_data["proof"] = {"edge": edge, "bound": 3-b}

            elif norm_C2 in masked_ref_vectors:
                 candidates = masked_ref_vectors[norm_C2]
                 match_data["identification"] = "implicit_subset_lower"
                 match_data["proof"] = {"subset": candidates[0]["subset"], "rhs": 3-b}

        matches.append(match_data)
        
        if match_data["identification"] == "unknown":
            active_edges = []
            for idx, val in enumerate(A):
                if val != 0:
                    active_edges.append((edges_ordered[idx], val))
            match_data["diagnostic_active_edges"] = str(active_edges)

    # Validate
    unknowns = [m for m in matches if m["identification"] == "unknown"]
    if unknowns:
        print(f"WARNING: {len(unknowns)} facets could not be identified.")
        for u in unknowns:
            print(f"  Facet {u['facet_id']} (b={u['ineq'][0]}): {u.get('diagnostic_active_edges')}")
    else:
        print("SUCCESS: All facets identified.")

    out_path = "RESULTS/facet_dictionary_n6_exact.json"
    with open(out_path, "w") as f:
        json.dump(matches, f, indent=2)
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    identify_facets()
