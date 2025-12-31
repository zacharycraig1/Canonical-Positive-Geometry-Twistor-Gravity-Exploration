import sys
import os
import json
from sage.all import *

# Add project root to path
sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents

def build_polytope_data():
    print("Phase G1: Building Forest Polytope Data...")
    
    # Ensure cache dir exists
    cache_dir = ".dcp_cache/phaseG"
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
        
    cases = [
        {"n": 4, "roots": [0, 1, 2]}, # Very small
        {"n": 5, "roots": [0, 1, 2]}, # Intermediate
        {"n": 6, "roots": [0, 1, 2]}  # The physical target
    ]
    
    summary = []
    
    for case in cases:
        n = case['n']
        roots = case['roots']
        label = f"n{n}_roots_{''.join(map(str, roots))}"
        print(f"\nProcessing {label}...")
        
        # 1. Get Exponents (Vertices)
        points, edges_ordered = get_forest_exponents(n, roots)
        print(f"  Found {len(points)} vertices.")
        
        # 2. Build Polytope
        P = Polyhedron(vertices=points)
        
        # 3. Compute Stats
        dim = P.dim()
        n_vertices = P.n_vertices()
        n_facets = P.n_facets()
        print(f"  Dim: {dim}, Vertices: {n_vertices}, Facets: {n_facets}")
        
        # 4. Extract Facet Inequalities
        # H-representation: A*x + b >= 0
        ieqs = P.inequalities()
        facets_data = []
        for ieq in ieqs:
            # ieq is [b, a1, a2, ...]
            # inequality is b + a.x >= 0
            # Convert sage Integers to python ints
            facets_data.append([int(x) for x in ieq])
            
        # 5. Save Data
        # Ensure vertices are also ints
        vertices_py = [[int(c) for c in v] for v in points]
        
        # Ensure edge tuples are ints
        edges_ordered_py = [[int(u), int(v)] for u, v in edges_ordered]

        data = {
            "n": int(n),
            "roots": [int(r) for r in roots],
            "dim": int(dim),
            "n_vertices": int(n_vertices),
            "n_facets": int(n_facets),
            "edges_ordered": edges_ordered_py,
            "vertices": vertices_py, 
            "facets_inequalities": facets_data
        }
        
        out_file = os.path.join(cache_dir, f"P_{label}.json")
        with open(out_file, "w") as f:
            json.dump(data, f, indent=2)
        print(f"  Saved to {out_file}")
        
        summary.append({
            "case": label,
            "vertices": n_vertices,
            "facets": n_facets,
            "dim": dim
        })
        
    print("\nPhase G1 Summary:")
    for s in summary:
        print(f"  {s['case']}: Dim {s['dim']}, V {s['vertices']}, F {s['facets']}")
        
if __name__ == "__main__":
    build_polytope_data()

