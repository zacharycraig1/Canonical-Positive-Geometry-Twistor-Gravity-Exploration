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

def generate_subset_vectors(n, edges_ordered):
    """
    Generates cut and internal vectors for all non-trivial subsets.
    Returns a dict: vector_tuple -> (type, subset)
    """
    vectors = {}
    
    # Subsets of size 2 to n-1 (since size 1 cuts are often trivial or same as size n-1)
    # Actually, iterate all masks from 1 to 2^n - 2
    # But for n=6, 2^6 = 64. Fast.
    
    edge_to_idx = {tuple(e): i for i, e in enumerate(edges_ordered)}
    num_edges = len(edges_ordered)
    
    for r in range(1, n): # size 1 to n-1
        for subset in itertools.combinations(range(n), r):
            S = set(subset)
            
            # 1. Cut Vector: edges between S and V\S
            cut_vec = [0] * num_edges
            for u, v in edges_ordered:
                idx = edge_to_idx[(u, v)]
                # Check if one in S and one not
                u_in = u in S
                v_in = v in S
                if u_in != v_in:
                    cut_vec[idx] = 1
            
            # Normalize and store
            norm_cut, _ = normalize_vector(cut_vec)
            if norm_cut not in vectors:
                vectors[norm_cut] = {"type": "cut", "subset": list(subset)}
                
            # 2. Internal Vector: edges within S
            if len(S) > 1:
                int_vec = [0] * num_edges
                for u, v in edges_ordered:
                    idx = edge_to_idx[(u, v)]
                    if u in S and v in S:
                        int_vec[idx] = 1
                
                norm_int, _ = normalize_vector(int_vec)
                if norm_int not in vectors:
                    vectors[norm_int] = {"type": "internal", "subset": list(subset)}

    return vectors

def match_facets():
    input_path = "RESULTS/facets_n6_ineq_exact.json"
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found. Run facet_audit_n6_strict.sage first.")
        return

    with open(input_path, "r") as f:
        data = json.load(f)
        
    n = data["n"]
    edges_ordered = data["edges_ordered"]
    facets = data["inequalities"] # List of [b, a0, a1, ...]
    
    print(f"Loaded {len(facets)} facets for n={n}")
    
    # Generate candidate vectors
    print("Generating subset vectors...")
    candidate_vectors = generate_subset_vectors(n, edges_ordered)
    print(f"Generated {len(candidate_vectors)} unique candidate vectors")
    
    # Match
    matches = []
    
    for i, facet_coeffs in enumerate(facets):
        # Parse coefficients
        # format: [b, A0, A1, ...]
        coeffs = [int(x) for x in facet_coeffs]
        b = coeffs[0]
        A = coeffs[1:]
        
        # Normalize A to compare
        norm_A, _ = normalize_vector(A)
        neg_norm_A, _ = normalize_vector([-x for x in A])
        
        match_info = None
        
        # Try direct match
        if norm_A in candidate_vectors:
            cand = candidate_vectors[norm_A]
            match_info = {
                "type": cand["type"],
                "subset": cand["subset"],
                "sign": +1,
                "aligned": True
            }
        elif neg_norm_A in candidate_vectors:
            cand = candidate_vectors[neg_norm_A]
            match_info = {
                "type": cand["type"],
                "subset": cand["subset"],
                "sign": -1,
                "aligned": True
            }
        else:
            match_info = {
                "type": "unknown",
                "subset": None,
                "aligned": False
            }
            
        entry = {
            "facet_id": i,
            "ineq_b_A": facet_coeffs,
            "best_match": match_info
        }
        matches.append(entry)
        
    # Stats
    matched_count = sum(1 for m in matches if m["best_match"]["aligned"])
    print(f"Matched {matched_count} / {len(facets)} facets to subsets.")
    
    output_path = "RESULTS/facet_dictionary_n6.json"
    with open(output_path, "w") as f:
        json.dump(matches, f, indent=2)
        
    print(f"Saved dictionary to {output_path}")

if __name__ == "__main__":
    match_facets()



