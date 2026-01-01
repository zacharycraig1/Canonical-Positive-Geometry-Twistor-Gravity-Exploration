import json
import os
import sys
# Ensure we can import from src
sys.path.append(os.getcwd())

from sage.all import Polyhedron, QQ
from src.posgeom.forest_polytope import get_forest_exponents

def audit_facets_n6():
    print("Generating forest exponents for n=6, roots=[0,1,2]...")
    # Use the existing function to get exponents (vertices) and edge ordering
    exponents, edges_ordered = get_forest_exponents(6, [0, 1, 2])
    
    print(f"Number of vertices (forests): {len(exponents)}")
    print(f"Dimension of ambient space: {len(edges_ordered)}")
    
    # Construct the polytope over QQ (exact arithmetic)
    print("Constructing Polytope over QQ...")
    P = Polyhedron(vertices=exponents, base_ring=QQ)
    
    print(f"Polytope dimension: {P.dim()}")
    print(f"Ambient dimension: {P.ambient_dim()}")
    
    # Get H-representation (inequalities)
    # The H-representation is a list of inequalities of the form: A*x + b >= 0
    # Sage returns them as HRepresentation objects
    print("Computing H-representation...")
    H = P.Hrepresentation()
    
    num_facets = len(H)
    print(f"Number of facets: {num_facets}")
    
    # Extract facets in a serializable format
    # Each inequality is usually stored as [b, a_0, a_1, ..., a_m] representing b + a*x >= 0
    # We want to store exactly what Sage gives us
    
    facets_list = []
    for h in H:
        # h.vector() returns (b, a_0, a_1, ...)
        # We ensure they are integers/rationals
        vec = list(h.vector())
        facets_list.append([str(x) for x in vec])
        
    results = {
        "n": int(6),
        "roots": [int(0), int(1), int(2)],
        "num_vertices": int(len(exponents)),
        "num_facets": int(num_facets),
        "edges_ordered": [[int(u), int(v)] for u, v in edges_ordered],
        "facets_b_A": facets_list,
        "note": "Each facet is [b, A_0, ..., A_m] where b + A.x >= 0"
    }
    
    # Ensure RESULTS directory exists
    os.makedirs("RESULTS", exist_ok=True)
    
    output_path = "RESULTS/facets_n6_exact.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
        
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    audit_facets_n6()

