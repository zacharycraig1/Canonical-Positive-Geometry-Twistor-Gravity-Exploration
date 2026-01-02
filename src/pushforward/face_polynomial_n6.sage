import json
import os
import sys
import itertools
from sage.all import *

# Ensure src imports work
sys.path.append(os.getcwd())

def get_spanning_trees_edges(nodes):
    """
    Return list of edge-sets for spanning trees of K_{nodes}.
    """
    k = len(nodes)
    if k == 1:
        return [[]]
    if k == 2:
        return [[tuple(sorted(tuple(nodes)))]]
    
    # Use Sage Graph for k >= 3
    # Construct complete graph on nodes
    G = Graph([nodes, list(itertools.combinations(nodes, 2))])
    trees = G.spanning_trees()
    
    forest_edge_sets = []
    for t in trees:
        edges = []
        for u, v, _ in t.edges():
            edges.append(tuple(sorted((u, v))))
        forest_edge_sets.append(edges)
    return forest_edge_sets

def enumerate_forests(n, roots):
    """
    Yield all rooted forests on n vertices with roots R.
    Each forest is a set of edges (tuples (u, v) with u < v).
    """
    non_roots = [r for r in range(n) if r not in roots]
    
    # Iterate over assignments of non-roots to roots
    for assignment in itertools.product(roots, repeat=len(non_roots)):
        comps = {r: [r] for r in roots}
        for i, owner in enumerate(assignment):
            comps[owner].append(non_roots[i])
            
        # For each component, we have a set of spanning trees.
        # Forest is a choice of one tree per component.
        
        comp_trees = []
        for r in roots:
            comp_trees.append(get_spanning_trees_edges(comps[r]))
            
        # Cartesian product of component trees
        for tree_tuple in itertools.product(*comp_trees):
            # tree_tuple is (edges_comp1, edges_comp2, edges_comp3)
            # Flatten to single edge set
            full_edges = []
            for t_edges in tree_tuple:
                full_edges.extend(t_edges)
            yield full_edges

def compute_face_polynomials():
    n = 6
    roots = [0, 1, 2]
    
    # Load boundary dictionary
    dict_path = os.path.join("RESULTS", f"boundary_dictionary_n{n}.json")
    if not os.path.exists(dict_path):
        print(f"Error: {dict_path} not found. Run R1 first.")
        return
        
    with open(dict_path, 'r') as f:
        boundary_dict = json.load(f)
        
    print("Enumerating forests...")
    all_forests = list(enumerate_forests(n, roots))
    print(f"Total forests: {len(all_forests)}")
    
    # Pre-compute edge indices for dot products
    # The inequalities use variable order from get_forest_polytope_inequalities
    # We need to match that order.
    # Reuse the function or replicate order logic:
    # (0,1), (0,2)... (0,5), (1,2)...
    edge_list = []
    for i in range(n):
        for j in range(i+1, n):
            edge_list.append((i, j))
    
    edge_to_idx = {e: k for k, e in enumerate(edge_list)}
    
    # Convert forests to incidence vectors (sparse or set)
    forest_sets = [set(f) for f in all_forests]
    
    # Prepare Polynomial Ring
    # z_ij vars
    var_names = [f"z_{u}_{v}" for u, v in edge_list]
    R = PolynomialRing(QQ, var_names)
    z_vars = R.gens()
    z_map = {edge_list[k]: z_vars[k] for k in range(len(edge_list))}
    
    results = {}
    
    for facet_id, info in boundary_dict.items():
        print(f"Processing facet {facet_id} ({info['type']})...")
        
        # Parse inequality
        # raw: [val, val, ...] corresponds to A.
        # We need to identify forests satisfying A.x == b?
        # Facets of P are defined by A.x <= b.
        # The facet ITSELF is the set where A.x == b.
        # So we check dot(A, incidence) == b. (using raw vals)
        # However, the stored dictionary has "inequality_raw" = A from A.x <= b ?
        # Let's check R1 script.
        # R1: "A_sage = vector(facet.A()); A_canon = -A_sage; b_canon = b_sage"
        # And "inequality_raw": [float(x) for x in A_canon]
        # So we have A_canon . x <= b_canon.
        # Facet is A_canon . x == b_canon.
        
        A = info["inequality_raw"]
        b = info["rhs"]
        
        # Find saturating forests
        saturating_forests = []
        for i, f_set in enumerate(forest_sets):
            # Compute dot product
            dot = 0
            for u, v in f_set:
                k = edge_to_idx[(u, v)]
                dot += A[k]
            
            # Check equality (floating point safe?)
            # A and b are from Sage QQ, but stored as float in JSON?
            # Better to be epsilon safe.
            if abs(dot - b) < 1e-9:
                saturating_forests.append(all_forests[i])
                
        print(f"  > Found {len(saturating_forests)} saturating forests.")
        
        # Build polynomial
        poly = R(0)
        for f_edges in saturating_forests:
            term = R(1)
            for e in f_edges:
                term *= z_map[e]
            poly += term
            
        # Factorize
        factored = poly.factor()
        
        # Store result
        results[facet_id] = {
            "type": info["type"],
            "num_forests": len(saturating_forests),
            "factorization": str(factored),
            "is_monomial": len(poly.monomials()) == 1 if poly != 0 else False,
            "is_zero": poly == 0
        }
        
    # Save results
    out_path = os.path.join("RESULTS", f"facet_factorizations_n{n}.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print(f"Saved factorizations to {out_path}")

if __name__ == "__main__":
    compute_face_polynomials()




