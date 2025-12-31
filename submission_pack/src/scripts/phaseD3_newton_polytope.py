import os
import sys
from sage.all import *

# Add src to path
sys.path.append(os.getcwd())

from src.chy_oracle.matrix_tree import get_k6_spanning_trees

def compute_newton_polytope():
    print("Computing Newton Polytope of the Spanning Tree Sum for K6...")
    
    trees = get_k6_spanning_trees()
    print(f"Number of trees: {len(trees)}")
    
    # Map edges (u,v) to indices 0..14
    edges = []
    for i in range(6):
        for j in range(i+1, 6):
            edges.append((i,j))
            
    edge_map = {e: k for k, e in enumerate(edges)}
    print(f"Number of edges (variables): {len(edges)}")
    
    # Create points
    points = []
    for tree in trees:
        pt = [0] * len(edges)
        for u, v in tree:
            if u > v: u, v = v, u
            idx = edge_map[(u,v)]
            pt[idx] = 1
        points.append(pt)
        
    print("Building Polyhedron (this may take a moment)...")
    P = Polyhedron(vertices=points)
    
    print(f"Dimension: {P.dim()}")
    print(f"Number of vertices: {P.n_vertices()}")
    print(f"Number of facets: {P.n_facets()}")
    print(f"Number of inequalities: {len(P.inequalities())}")
    
    # Check if it matches Spanning Tree Polytope properties
    # Dimension should be |E| - |V| + 1 ? No.
    # Spanning tree polytope lies in the hyperplane sum x_e = n-1.
    # Dim = |E| - 1. (Actually depends on correlations).
    # For K6, |E|=15. |V|=6. Tree size = 5.
    # Vertices are in 15D space.
    # Constraints: x_e >= 0. Sum x_e = 5.
    # Plus subtour elimination constraints.
    
    print("\nSample Facet Inequalities:")
    for i, ineq in enumerate(P.inequalities()[:5]):
        print(f"  {ineq}")
        
    # Check if Permutohedron?
    # Spanning tree polytope of Kn is a generalized permutohedron.
    
if __name__ == "__main__":
    compute_newton_polytope()

