import sys
import os
from sage.all import *

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.forest_sum import k_rooted_forests_complete_graph

def analyze_polytopes():
    print("Starting Phase F3.1: Newton Polytope Analysis...")
    
    n = 6
    
    # ---------------------------------------------------------
    # 1. Spanning Tree Polytope (Sanity Check)
    # ---------------------------------------------------------
    print(f"\n[1/2] Computing Spanning Tree Polytope for K_{n}...")
    K6 = graphs.CompleteGraph(n)
    trees = list(K6.spanning_trees())
    print(f"  Found {len(trees)} spanning trees.")
    
    # Edges canonical order
    edges_ordered = []
    edge_to_idx = {}
    counter = 0
    for i in range(n):
        for j in range(i + 1, n):
            edges_ordered.append((i, j))
            edge_to_idx[(i, j)] = counter
            counter += 1
    dim_ambient = len(edges_ordered)
    
    tree_points = []
    for tree in trees:
        pt = [0] * dim_ambient
        for u, v, _ in tree.edges():
            if u > v: u, v = v, u
            idx = edge_to_idx[(u, v)]
            pt[idx] = 1
        tree_points.append(pt)
        
    P_tree = Polyhedron(vertices=tree_points)
    print(f"  Tree Polytope: Dim {P_tree.dim()}, Vertices {P_tree.n_vertices()}, Facets {P_tree.n_facets()}")
    
    # Standard Permutohedron check: Dim should be n-1? No, embedded in edge space.
    # Spanning Tree Polytope of K_n has dimension n(n-1)/2 - n? 
    # Actually it is a projection of permutohedron.
    # Dim is n-1?
    # Let's save stats.
    
    # ---------------------------------------------------------
    # 2. Rooted Forest Polytope (The Main Object)
    # ---------------------------------------------------------
    roots = [0, 1, 2]
    print(f"\n[2/2] Computing 3-Rooted Forest Polytope for K_{n} (roots: {roots})...")
    
    forests = list(k_rooted_forests_complete_graph(n, roots))
    print(f"  Found {len(forests)} forests.")
    
    forest_points = []
    for forest in forests:
        pt = [0] * dim_ambient
        for u, v in forest:
            if u > v: u, v = v, u
            idx = edge_to_idx[(u, v)]
            pt[idx] = 1
        forest_points.append(pt)
        
    P_forest = Polyhedron(vertices=forest_points)
    print(f"  Forest Polytope: Dim {P_forest.dim()}, Vertices {P_forest.n_vertices()}, Facets {P_forest.n_facets()}")
    
    # Save results
    results = {
        "n": n,
        "tree_polytope": {
            "vertices": int(P_tree.n_vertices()),
            "facets": int(P_tree.n_facets()),
            "dim": int(P_tree.dim())
        },
        "forest_polytope": {
            "vertices": int(P_forest.n_vertices()),
            "facets": int(P_forest.n_facets()),
            "dim": int(P_forest.dim())
        }
    }
    
    import json
    with open("results/phaseF_newton_polytopes.json", "w") as f:
        json.dump(results, f, indent=2)
        
    save(P_forest, "results/forest_polytope_n6_k3.sobj")
    
    print("\nAnalysis Complete. JSON saved.")

if __name__ == "__main__":
    analyze_polytopes()
