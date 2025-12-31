import sys
import os
from sage.all import *
import itertools

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian

def k_rooted_forests_complete_graph(n, roots):
    """
    Yields spanning forests of K_n with exactly |roots| components,
    where each component contains exactly one root from 'roots'.
    
    Each forest is returned as a list of edges [(u,v), ...].
    
    Algorithm: Brute force edge combinations.
    Number of edges must be n - |roots|.
    """
    num_roots = len(roots)
    num_edges = n - num_roots
    
    # All edges in K_n
    all_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            all_edges.append((i, j))
            
    # Iterate over combinations
    for edges in itertools.combinations(all_edges, num_edges):
        # Check 1: No cycles (must be a forest)
        # Check 2: Connectivity to roots
        
        # Build adjacency for this set of edges
        adj = {i: [] for i in range(n)}
        for u, v in edges:
            adj[u].append(v)
            adj[v].append(u)
            
        # Find components
        visited = set()
        components = []
        is_forest = True
        
        # We can check cycle count: V - E = K. 
        # Here n - (n-k) = k components if acyclic.
        # So we just need to check if we have exactly k components 
        # and each component has exactly 1 root.
        
        for i in range(n):
            if i not in visited:
                comp_roots = 0
                stack = [i]
                visited.add(i)
                comp_nodes = []
                
                while stack:
                    curr = stack.pop()
                    comp_nodes.append(curr)
                    if curr in roots:
                        comp_roots += 1
                    
                    for neighbor in adj[curr]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)
                        # Note: DFS cycle check requires parent pointer, 
                        # but simple component counting is sufficient here
                        # because |E| = n - k. If components = k, it must be a forest.
                
                components.append(comp_nodes)
                
                # Optimization: Fail early if a component has 0 or >1 roots
                if comp_roots != 1:
                    is_forest = False
                    break
        
        if is_forest and len(components) == num_roots:
            yield edges

def forest_sum_minor(n, lambdas, tildes, x, y, roots=(0, 1, 2)):
    """
    Computes the sum over rooted spanning forests of the product of edge weights.
    Weight for edge (i,j) is a_{ij} = w_{ij} * C_i * C_j.
    
    Corresponds to det(L_tilde^(roots)).
    """
    # Precompute edge weights
    # w_ij = [ij]/<ij>
    # C_i = <ix><iy>
    
    weights = {}
    Cs = []
    
    for i in range(n):
        c_val = ang_bracket(lambdas[i], x) * ang_bracket(lambdas[i], y)
        Cs.append(c_val)
        
    for i in range(n):
        for j in range(i+1, n):
            ang = ang_bracket(lambdas[i], lambdas[j])
            sq = sq_bracket(tildes[i], tildes[j])
            if ang == 0:
                raise ValueError(f"Angle bracket <{i}{j}> is zero")
            w_ij = sq / ang
            
            # The Laplacian has -w_ij C_i C_j on off-diagonal.
            # Matrix Tree Theorem sums product of (-L_ij) = product of (w_ij C_i C_j).
            weights[(i, j)] = w_ij * Cs[i] * Cs[j]
            
    total_sum = 0
    
    # Iterate forests
    for forest in k_rooted_forests_complete_graph(n, roots):
        term = 1
        for u, v in forest:
            # Order u < v for lookup
            if u > v: u, v = v, u
            term *= weights[(u, v)]
        total_sum += term
        
    return total_sum

def verify_forest_theorem_n6():
    print("Verifying Forest Sum Theorem for n=6...")
    
    # Sample kinematics
    from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
    lambdas, tildes = sample_spinors_from_twistor(seed=0, n=6)
    
    # Reference spinors
    x = vector(QQ, [1, 2])
    y = vector(QQ, [3, 1])
    
    # 1. Compute Determinant
    L_tilde, _, _ = hodges_weighted_laplacian(lambdas, tildes, x, y)
    indices = [3, 4, 5] # Rows to keep (roots are 0,1,2, so delete them)
    # Wait, roots are the deleted rows/cols in the minor
    det_val = L_tilde.matrix_from_rows_and_columns(indices, indices).det()
    
    # 2. Compute Forest Sum
    forest_val = forest_sum_minor(6, lambdas, tildes, x, y, roots=[0, 1, 2])
    
    print(f"Determinant: {det_val}")
    print(f"Forest Sum:  {forest_val}")
    
    if abs(det_val - forest_val) < 1e-10:
        print("SUCCESS: Exact match.")
    else:
        print("FAIL: Mismatch.")
        
if __name__ == "__main__":
    verify_forest_theorem_n6()

