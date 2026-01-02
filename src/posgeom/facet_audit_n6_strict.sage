import json
import os
import sys

# Ensure we can import from src
sys.path.append(os.getcwd())

from sage.all import Polyhedron, QQ
from src.posgeom.forest_polytope import get_forest_exponents

def audit_facets_n6_strict():
    print("Generating forest exponents for n=6, roots=[0,1,2]...")
    exponents, edges_ordered = get_forest_exponents(6, [0, 1, 2])
    
    print(f"Number of vertices (forests): {len(exponents)}")
    print(f"Dimension of ambient space: {len(edges_ordered)}")
    
    # Construct the polytope over QQ (exact arithmetic)
    print("Constructing Polytope over QQ...")
    P = Polyhedron(vertices=exponents, base_ring=QQ)
    
    print(f"Polytope dimension: {P.dim()}")
    print(f"Ambient dimension: {P.ambient_dim()}")
    
    # Strict separation
    inequalities = P.inequalities()
    equations = P.equations()
    
    print(f"Number of inequalities (facets): {len(inequalities)}")
    print(f"Number of equations: {len(equations)}")
    
    # Helper to serialize H-representation
    def serialize_hrep(h_list):
        serialized = []
        for h in h_list:
            # h.vector() is [b, a0, a1, ...] for b + A.x >= 0 or == 0
            vec = list(h.vector())
            serialized.append([str(x) for x in vec])
        return serialized
    
    facets_data = {
        "n": 6,
        "roots": [0, 1, 2],
        "dim": int(P.dim()),
        "ambient_dim": int(P.ambient_dim()),
        "count": len(inequalities),
        "inequalities": serialize_hrep(inequalities),
        "edges_ordered": [[int(u), int(v)] for u, v in edges_ordered],
        "note": "Each entry is [b, A_0, ..., A_m] where b + A.x >= 0"
    }
    
    eq_data = {
        "n": 6,
        "roots": [0, 1, 2],
        "dim": int(P.dim()),
        "ambient_dim": int(P.ambient_dim()),
        "count": len(equations),
        "equations": serialize_hrep(equations),
        "edges_ordered": [[int(u), int(v)] for u, v in edges_ordered],
        "note": "Each entry is [b, A_0, ..., A_m] where b + A.x == 0"
    }
    
    os.makedirs("RESULTS", exist_ok=True)
    
    with open("RESULTS/facets_n6_ineq_exact.json", "w") as f:
        json.dump(facets_data, f, indent=2)
    print("Saved RESULTS/facets_n6_ineq_exact.json")
        
    with open("RESULTS/facets_n6_eq_exact.json", "w") as f:
        json.dump(eq_data, f, indent=2)
    print("Saved RESULTS/facets_n6_eq_exact.json")

if __name__ == "__main__":
    audit_facets_n6_strict()




