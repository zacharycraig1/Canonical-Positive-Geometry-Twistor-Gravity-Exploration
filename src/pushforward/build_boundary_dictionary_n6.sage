import json
import os
import sys
# Ensure we can import from src
sys.path.append(os.getcwd())
from sage.all import *
from src.posgeom.forest_polytope_inequalities import get_forest_polytope_inequalities

def classify_facet(ineq, edges):
    """
    Classify a facet inequality: A*x <= b
    ineq: Sage inequality object or (A, b) tuple
    edges: list of (u, v) corresponding to indices of A
    
    Returns: dict with classification info
    """
    if hasattr(ineq, 'A'):
        A = ineq.A()
        b = ineq.b()
    else:
        # Tuple (vector, rhs)
        A = ineq[0]
        b = ineq[1]
    
    # Analyze A
    # Check for unit vectors (bounds)
    non_zeros = [(i, val) for i, val in enumerate(A) if val != 0]
    
    classification = {
        "type": "unknown",
        "rhs": float(b),
        "inequality_raw": [float(x) for x in A]
    }
    
    if len(non_zeros) == 1:
        idx, val = non_zeros[0]
        edge = edges[idx]
        classification["edge"] = [int(edge[0]), int(edge[1])]
        
        # x_e >= 0  <=> -x_e <= 0  (val = -1, b = 0)
        if val == -1 and b == 0:
            classification["type"] = "lower_bound"
            classification["variable"] = f"x_{edge[0]}_{edge[1]}"
            classification["meaning"] = "edge_vanishes"
            
        # x_e <= 1  <=> x_e <= 1 (val = 1, b = 1)
        elif val == 1 and b == 1:
            classification["type"] = "upper_bound"
            classification["variable"] = f"x_{edge[0]}_{edge[1]}"
            classification["meaning"] = "edge_saturated"
    else:
        # Check for subset rank type
        # All coeffs should be 0 or 1 (or close to)
        is_subset_sum = all(val == 1 for _, val in non_zeros)
        if is_subset_sum:
            # Reconstruct subset S
            # The edges with coeff 1 are exactly E(S)
            involved_nodes = set()
            involved_edges = []
            for idx, _ in non_zeros:
                u, v = edges[idx]
                involved_nodes.add(u)
                involved_nodes.add(v)
                involved_edges.append((u, v))
            
            # Verify if involved_edges == E(involved_nodes)
            # Actually, for the inequality x(E(S)) <= rhs, the edges MUST be exactly all edges in S?
            # Yes, if it's a standard subset inequality.
            
            S = sorted(list(involved_nodes))
            classification["type"] = "subset_rank"
            classification["subset_S"] = [int(x) for x in S]
            classification["roots_in_S"] = len(set(S).intersection({0, 1, 2})) # Assuming roots 0,1,2
            classification["num_edges"] = len(involved_edges)
            classification["meaning"] = "subset_factorization"
            
    return classification

def build_boundary_dictionary(n=6, roots=[0, 1, 2]):
    print(f"Generating inequalities for n={n}, roots={roots}...")
    edges, ineqs_list = get_forest_polytope_inequalities(n, roots)
    
    print(f"Building Polyhedron from {len(ineqs_list)} inequalities...")
    # Sage Polyhedron expects ieqs as list of [b, -A] for A*x + b >= 0 ?? 
    # Or [b, a1, a2...] for b + A*x >= 0 ?
    # Let's check Sage docs or use interactive definition.
    # Usually: ieqs=[[b, a1, ..., an], ...] for b + a.x >= 0
    # Our format: A.x <= b  =>  b - A.x >= 0
    # So [b, -A]
    
    sage_ieqs = []
    for vec, rhs in ineqs_list:
        # [rhs, -vec[0], -vec[1], ...]
        line = [rhs] + [-x for x in vec]
        sage_ieqs.append(line)
        
    P = Polyhedron(ieqs=sage_ieqs)
    print(f"Polyhedron built. Vertices: {P.n_vertices()}, Facets: {P.n_facets()}")
    
    # We want to map FACETS.
    # P.inequalities() gives the minimal set of inequalities defining facets.
    facets = P.inequalities()
    
    dictionary = {}
    
    for i, facet in enumerate(facets):
        # facet is an Inequality object
        # It corresponds to b + A.x >= 0
        # We want to convert back to A.x <= b form to match our classification logic
        # A_sage = facet.A()
        # b_sage = facet.b()
        # b + A.x >= 0  <=>  -A.x <= b
        # So our A_canonical = -A_sage, b_canonical = b_sage
        
        A_sage = vector(facet.A())
        b_sage = facet.b()
        
        A_canon = -A_sage
        b_canon = b_sage
        
        # We classify based on A_canon * x <= b_canon
        cls = classify_facet((A_canon, b_canon), edges)
        dictionary[i] = cls
        
    # Output to JSON
    output_path = os.path.join("RESULTS", f"boundary_dictionary_n{n}.json")
    os.makedirs("RESULTS", exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(dictionary, f, indent=2)
        
    print(f"Boundary dictionary saved to {output_path}")
    return dictionary

if __name__ == "__main__":
    build_boundary_dictionary()

