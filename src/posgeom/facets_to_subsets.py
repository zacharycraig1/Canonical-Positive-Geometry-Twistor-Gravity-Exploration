import json
import itertools
import sys
import os

# Ensure src is in path so we can import from posgeom
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from sage.all import *
from posgeom.forest_polytope import get_forest_exponents

def generate_facets_n6():
    n = 6
    roots = [0, 1, 2]
    
    print(f"Generating forest exponents for n={n}, roots={roots}...")
    exponents, edge_order = get_forest_exponents(n, roots)
    print(f"Number of forests (vertices): {len(exponents)}")
    
    print("Computing Convex Hull (this might take a few seconds)...")
    P = Polyhedron(vertices=exponents, base_ring=QQ)
    
    facets = P.Hrepresentation()
    print(f"Number of facets: {len(facets)}")
    
    return facets, edge_order

def classify_facet(ineq, constant, edge_order, n, roots):
    """
    Classifies a facet inequality: sum(c_e * x_e) + constant >= 0.
    
    Returns a dictionary description.
    """
    # Inequality vector 'ineq' corresponds to edges in 'edge_order'
    
    # 1. Check for trivial bounds: x_e >= 0 or x_e <= 1
    # x_e >= 0  =>  1*x_e + 0 >= 0
    # x_e <= 1  => -1*x_e + 1 >= 0
    
    non_zeros = [(i, c) for i, c in enumerate(ineq) if c != 0]
    
    if len(non_zeros) == 1:
        idx, coeff = non_zeros[0]
        edge = edge_order[idx]
        if coeff == 1 and constant == 0:
            return {
                "type": "lower_bound",
                "edge": edge,
                "desc": f"x_{edge} >= 0"
            }
        elif coeff == -1 and constant == 1:
            return {
                "type": "upper_bound",
                "edge": edge,
                "desc": f"x_{edge} <= 1"
            }
            
    # 2. Check for Subset Inequalities (Upper Bounds)
    # Form: sum_{e in E(S)} x_e <= |S| - |R cap S|
    # rewritten: -sum x_e + (|S| - |R cap S|) >= 0
    # Coeffs are -1.
    
    is_upper = True
    for idx, c in non_zeros:
        if c != -1:
            is_upper = False
            break
            
    if is_upper and constant > 0:
        involved_edges = [edge_order[idx] for idx, c in non_zeros]
        S_vertices = set()
        for u, v in involved_edges:
            S_vertices.add(u)
            S_vertices.add(v)
            
        # Add roots to S if it makes the set of edges match?
        # Actually, the subset inequality usually involves a set S.
        # The edges are E(S).
        # Since root-root edges are 0, they won't appear in involved_edges.
        # But if S contains roots, we should expect ALL edges between S-nodes to be present (except root-root).
        
        # We need to guess S.
        # Hypothesis: S is exactly the set of vertices incident to the involved edges.
        # OR: S might include some roots that are isolated in the involved graph but necessary for the formula?
        # No, if a vertex is in S, and we sum over E(S), we expect edges connecting it to be involved.
        
        s_list = sorted(list(S_vertices))
        
        # Expected edges in E(S) excluding root-root edges
        expected_edges = []
        for i in range(len(s_list)):
            for j in range(i+1, len(s_list)):
                u, v = s_list[i], s_list[j]
                if u > v: u, v = v, u
                
                # Exclude root-root edges
                if u in roots and v in roots:
                    continue
                    
                expected_edges.append((u, v))
                
        # Check match
        if set(involved_edges) == set(expected_edges):
            r_cap_s = len([r for r in roots if r in S_vertices])
            # If S has roots, we need |S|-|R_in_S| edges max (one comp per root).
            # If S has NO roots, we need |S|-1 edges max (no cycles).
            rhs = len(S_vertices) - max(1, r_cap_s)
            
            if constant == rhs:
                return {
                    "type": "subset_sum_upper",
                    "subset": list(s_list),
                    "rhs": rhs,
                    "desc": f"sum(x_e for e in S) <= {rhs}, S={list(s_list)}"
                }
            elif constant == (n - len(roots)) and set(involved_edges) == set(edge_order):
                 return {
                    "type": "global_upper",
                    "desc": f"sum(all x_e) <= {constant}"
                 }
            else:
                 return {
                    "type": "subset_mismatch",
                    "subset": list(s_list),
                    "expected_rhs": rhs,
                    "actual_rhs": constant,
                    "desc": f"Subset S={list(s_list)} match, but RHS {constant} != {rhs}"
                }
        
        # Check for Implicit Non-Negativity (sum_{others} <= GlobalMax)
        # This occurs if x_miss >= 0 is a facet, but x_miss = GlobalMax - sum_{others}.
        # So sum_{others} <= GlobalMax.
        global_max = n - len(roots)
        if constant == global_max:
             # Check if involved_edges is exactly ALL edges minus ONE edge
             all_edges_set = set(edge_order)
             inv_set = set(involved_edges)
             diff = all_edges_set - inv_set
             
             # Filter out root-root edges from diff
             real_diff = [e for e in diff if not (e[0] in roots and e[1] in roots)]
             
             if len(real_diff) == 1:
                 missing = real_diff[0]
                 return {
                     "type": "implicit_non_negativity",
                     "missing_edge": str(missing),
                     "desc": f"x_{missing} >= 0 (implicit via global sum)"
                 }

    # 3. Check for Lower Bounds (Positive Coeffs)
    # Form: sum c_e x_e - k >= 0  => sum x_e >= k
    is_lower = True
    for idx, c in non_zeros:
        if c != 1:
            is_lower = False
            break
            
    if is_lower and constant < 0:
        k = -constant
        involved_edges = [edge_order[idx] for idx, c in non_zeros]
        
        # Check if it's "sum incident to v >= 1"
        # Find common vertex?
        # Or Global Sum >= 3
        if set(involved_edges) == set(edge_order) and k == (n - len(roots)):
             return {
                "type": "global_lower",
                "desc": f"sum(all x_e) >= {k}"
             }
        
        # Check degree constraint
        # Degree of v >= 1. Involved edges are exactly those incident to v.
        # And v is non-root.
        counts = {}
        for u,v in involved_edges:
            counts[u] = counts.get(u,0)+1
            counts[v] = counts.get(v,0)+1
            
        # If it is degree constraint, there is a central vertex v connected to all other u in involved_edges.
        # And involved_edges contains ALL valid edges incident to v.
        
        # Let's just dump it as lower bound sum
        return {
            "type": "sum_lower",
            "k": k,
            "edges": [str(e) for e in involved_edges],
            "desc": f"sum(subset) >= {k}"
        }

    # If falling through
    return {
        "type": "unknown",
        "coeffs": {str(edge_order[i]): str(c) for i, c in non_zeros},
        "constant": str(constant),
        "desc": "Unclassified"
    }

def sage_default(obj):
    if isinstance(obj, (Integer, int)):
        return int(obj)
    if isinstance(obj, (Rational, float)):
        return float(obj)
    return str(obj)

def main():
    facets, edge_order = generate_facets_n6()
    n = 6
    roots = [0, 1, 2]
    
    results = []
    
    print("\nClassifying Facets...")
    for f in facets:
        # Sage H-rep: A*x + b >= 0
        # f.A() is vector of coeffs, f.b() is constant
        # Note: Sage might return sparse vector.
        ineq = list(f.A())
        constant = f.b()
        
        info = classify_facet(ineq, constant, edge_order, n, roots)
        results.append(info)
        print(f"  {info['desc']}")
        if info['type'] == 'unknown':
             print(f"    Coeffs: {info['coeffs']}, Const: {info['constant']}")
        
    # Stats
    counts = {}
    for r in results:
        t = r['type']
        counts[t] = counts.get(t, 0) + 1
        
    print("\nSummary:")
    for t, c in counts.items():
        print(f"  {t}: {c}")
        
    # Save to JSON
    out_path = "phase_u/facet_dictionary_n6.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=sage_default)
    print(f"\nSaved dictionary to {out_path}")

if __name__ == "__main__":
    main()

