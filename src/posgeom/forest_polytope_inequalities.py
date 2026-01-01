import sys
from sage.all import *
from itertools import combinations

def get_forest_polytope_inequalities(n, roots):
    """
    Generate inequalities for the Forest Polytope (spanning forests of Kn rooted at 'roots').
    This is a face of the Spanning Tree Polytope of Kn + super-root.
    
    Variables: x_{ij} for 0 <= i < j < n.
    Total variables: n*(n-1)/2.
    
    Constraints:
    1. 0 <= x_{ij} <= 1
    2. Sum of all edges = n - |roots| (since |roots| are connected to super-root)
    3. For any subset S of vertices, sum(x_e for e in S) <= |S| - 1
       (Wait, if S contains ALL roots, it can be larger?
        No, the condition for forests is:
        For any subset S, x(E(S)) <= |S| - k_S where k_S is number of connected components?
        Or simply: x(E(S)) <= |S| - 1 is for trees.
        For forests rooted at R:
        If S \cap R is empty: x(E(S)) <= |S| - 1.
        If S \cap R is not empty: x(E(S)) <= |S| - |S \cap R|? No.
        
        Let's use the Super-Root construction.
        G' = Kn + {rho}. Edges (rho, r) for r in roots.
        x_{rho, r} = 1 for r in roots.
        x_{rho, v} = 0 for v not in roots.
        
        Tree Polytope constraints on G':
        Sum_{e in G'} x_e = (n+1) - 1 = n.
        For any S' subset of V(G'), x(E(S')) <= |S'| - 1.
        
        Substitute x_{rho, r} = 1:
        Total sum on Kn edges: n - |roots|.
        
        Subsets S subset V(Kn):
        Case 1: S does not involve rho.
           x(E(S)) <= |S| - 1.
        Case 2: S' = S + {rho}.
           x(E(S')) <= |S| + 1 - 1 = |S|.
           E(S') includes E(S) and edges (rho, v) for v in S.
           So x(E(S)) + sum_{v in S \cap roots} 1 <= |S|.
           x(E(S)) <= |S| - |S \cap roots|.
           
        So for any S subset V(Kn):
           x(E(S)) <= |S| - (1 if S \cap roots is empty else |S \cap roots| ?)
           Wait.
           If S \cap roots is empty: x(E(S)) <= |S| - 1. (Standard forest cond)
           If S \cap roots is not empty:
             x(E(S)) <= |S| - |S \cap roots|.
             
        Example: S = {r1, r2} (two roots).
        x(r1, r2) <= 2 - 2 = 0.
        So no edges between roots! That's correct for rooted forests (paths go to roots, not between them?).
        Actually, in matrix-tree theorem, roots are where edges go?
        If we have edges between roots, we might form a cycle if we view them as connected at super-root.
        Yes, in G', roots are connected to rho. If r1-r2 exists, r1-rho-r2 is a cycle.
        So no edges between roots. Correct.
    
    Output: List of (A, b) where A.x <= b.
    """
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))
            
    num_vars = len(edges)
    inequalities = []
    
    # x_e >= 0 -> -x_e <= 0
    for k in range(num_vars):
        vec = vector(QQ, num_vars)
        vec[k] = -1
        inequalities.append((vec, 0))
        
    # x_e <= 1 -> x_e <= 1
    for k in range(num_vars):
        vec = vector(QQ, num_vars)
        vec[k] = 1
        inequalities.append((vec, 1))
        
    # Subsets S
    # Iterate all 2^n - 1 subsets. For n=6, 63 subsets. Fast.
    for r in range(1, n+1):
        for S in combinations(range(n), r):
            S_set = set(S)
            
            # Count roots in S
            roots_in_S = len(S_set.intersection(set(roots)))
            
            rhs = 0
            if roots_in_S == 0:
                rhs = len(S) - 1
            else:
                rhs = len(S) - roots_in_S
                
            # Build LHS vector (sum of x_e for e in E(S))
            vec = vector(QQ, num_vars)
            
            # Edges in S
            for k, (u, v) in enumerate(edges):
                if u in S_set and v in S_set:
                    vec[k] = 1
                    
            # Check if this inequality is trivial (e.g. 0 <= 0 or if no edges)
            if vec.is_zero():
                if rhs < 0:
                    print(f"Error: Impossible constraint for S={S}, rhs={rhs}")
                continue
                
            inequalities.append((vec, rhs))
            
    # Equality constraint: Sum x_e = n - |roots|
    # WAIT. The super-root construction for rooted forests:
    # We consider forests F on Kn such that each component contains exactly one root.
    # The number of edges in such a forest is n - |roots|.
    # But wait. n vertices. |roots| components.
    # Total vertices = n. Components = k. Edges = n - k.
    # Correct.
    
    # Are there other equalities?
    # In the Spanning Tree Polytope of G', sum x_e = n.
    # x_{rho, r} = 1 fixed.
    # So sum_{e in Kn} x_e = n - |roots|.
    
    # But maybe we are setting x_{rho, v} = 0 for v not in roots implies 
    # we removed those edges.
    
    # Let's check why dim dropped to 11 (from 14).
    # 15 vars. 1 equality -> 14.
    # We got 11. Lost 3 dimensions.
    # 3 roots? Maybe x_{r1, r2} = 0 for roots?
    # My manual logic earlier: "So no edges between roots!"
    # x(r1, r2) <= 0.
    # Since x >= 0, this implies x(r1, r2) = 0.
    # If roots are {0, 1, 2}, edges (0,1), (0,2), (1,2) are forced to 0.
    # 3 edges forced to 0.
    # 14 - 3 = 11.
    # Matches exactly!
    # So the polytope lives in a lower dimensional face where edges between roots are 0.
    
    total_rhs = n - len(roots)
    total_vec = vector(QQ, [1]*num_vars)
    
    inequalities.append((total_vec, total_rhs))
    inequalities.append((-total_vec, -total_rhs))
    
    return edges, inequalities

if __name__ == "__main__":
    import argparse
    
    # Parse args manually if needed or just run for n=6
    n = 6
    roots = [0, 1, 2]
    
    print(f"Generating inequalities for n={n}, roots={roots}...")
    edges, ineqs = get_forest_polytope_inequalities(n, roots)
    
    print(f"Number of variables: {len(edges)}")
    print(f"Number of inequalities: {len(ineqs)}")
    
    # Save to JSON or similar if needed, or just print sample
    # We will use this module in compare_facets.py

