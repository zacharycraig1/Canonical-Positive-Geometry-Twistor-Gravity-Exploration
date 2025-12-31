import sys
import os
from sage.all import *

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.forest_sum import k_rooted_forests_complete_graph

def analyze_forest_polytope():
    print("Starting Phase F2b: Forest Newton Polytope Analysis...")
    
    n = 6
    roots = [0, 1, 2]
    print(f"Generating 3-rooted forests for K_{n} (roots: {roots})...")
    
    # 1. Enumerate Forests
    forests = list(k_rooted_forests_complete_graph(n, roots))
    num_forests = len(forests)
    print(f"Found {num_forests} forests.")
    
    # 2. Map edges to coordinates
    # We need a canonical ordering of edges in K_n
    edges_ordered = []
    edge_to_idx = {}
    counter = 0
    for i in range(n):
        for j in range(i + 1, n):
            edges_ordered.append((i, j))
            edge_to_idx[(i, j)] = counter
            counter += 1
            
    dim_ambient = len(edges_ordered)
    print(f"Ambient dimension (edges): {dim_ambient}")
    
    # 3. Create vertices
    points = []
    for forest in forests:
        pt = [0] * dim_ambient
        for u, v in forest:
            if u > v: u, v = v, u
            idx = edge_to_idx[(u, v)]
            pt[idx] = 1
        points.append(pt)
        
    # 4. Build Polytope
    print("Constructing Polytope (this may take a moment)...")
    P = Polyhedron(vertices=points)
    
    print(f"Polytope Dimension: {P.dim()}")
    print(f"Vertices: {P.n_vertices()}")
    print(f"Facets: {P.n_facets()}")
    
    # 5. Save Object
    save_path = "forest_polytope_n6_k3.sobj"
    save(P, save_path)
    print(f"Saved polytope to {save_path}")
    
    # 6. Comparison with Spanning Tree Polytope
    # Spanning trees have n-1 edges. These have n-3 edges.
    # This is a generalized permutohedron.
    
    # Check if it's a "simple" polytope?
    # P.is_simple() # might be slow
    
    print("\nAnalysis Complete.")

if __name__ == "__main__":
    analyze_forest_polytope()

