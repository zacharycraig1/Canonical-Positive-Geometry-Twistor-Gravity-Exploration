import os
import sys
from sage.all import *

# Fix path for imports when run as script
if __name__ == "__main__":
    sys.path.append(os.getcwd())

try:
    from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket, sq_bracket
except ImportError:
    from chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket, sq_bracket

# Cache for spanning trees of K6
_K6_SPANNING_TREES = None

def get_k6_spanning_trees():
    """
    Returns a list of spanning trees for the complete graph K6.
    Each tree is represented as a list of edges (tuples (i, j)).
    """
    global _K6_SPANNING_TREES
    if _K6_SPANNING_TREES is not None:
        return _K6_SPANNING_TREES
    
    # Generate K6
    K6 = graphs.CompleteGraph(6)
    # Get spanning trees
    trees = K6.spanning_trees()
    
    # Convert to list of edges for faster iteration
    _K6_SPANNING_TREES = [tree.edges(labels=False) for tree in trees]
    return _K6_SPANNING_TREES

def laplacian_from_weights(weights, n):
    """
    Constructs the n x n Laplacian matrix L from a dictionary of edge weights.
    L_ij = -w_ij  for i != j
    L_ii = sum_{k!=i} w_ik
    
    Args:
        weights: dict {(i,j): weight} for i < j
        n: number of vertices
    """
    L = matrix(QQ, n, n)
    
    # Helper to get weight safely
    def get_w(i, j):
        if i == j: return 0
        if i > j: return weights.get((j, i), 0)
        return weights.get((i, j), 0)
        
    for i in range(n):
        row_sum = 0
        for j in range(n):
            if i != j:
                w = get_w(i, j)
                L[i, j] = -w
                row_sum += w
        L[i, i] = row_sum
        
    return L

def tree_sum_kirchhoff(weights, n, delete=0):
    """
    Computes the weighted sum of spanning trees using the Matrix-Tree Theorem (Kirchhoff).
    Returns det(L^(delete)), the principal minor removing row/col 'delete'.
    """
    L = laplacian_from_weights(weights, n)
    
    # Remove row and column 'delete'
    indices = [i for i in range(n) if i != delete]
    L_minor = L.matrix_from_rows_and_columns(indices, indices)
    
    return L_minor.det()

def hodges_weighted_laplacian(lambdas, tilde_lambdas, x, y):
    """
    Constructs the weighted Laplacian corresponding to Hodges' formula.
    
    Weights: w_ij = [ij]/<ij>
    Vertex weights: C_i = <i x><i y>
    
    Weighted Laplacian L_tilde:
      L_tilde_ij = -w_ij * C_i * C_j
      L_tilde_ii = sum_{k!=i} w_ik * C_i * C_k
      
    Returns:
        L_tilde: The weighted Laplacian matrix
        C: List of vertex weights [C_0, ..., C_{n-1}]
        Phi_tilde: The matrix -D^{-1} L_tilde D^{-1} which resembles Hodges Phi
    """
    n = len(lambdas)
    
    # 1. Compute basic weights w_ij
    weights = {}
    for i in range(n):
        for j in range(i + 1, n):
            ang = ang_bracket(lambdas[i], lambdas[j])
            sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
            if ang == 0:
                raise ValueError(f"Angle bracket <{i}{j}> is zero")
            weights[(i, j)] = sq / ang
            
    # 2. Compute vertex weights C_i
    C = []
    for i in range(n):
        c_val = ang_bracket(lambdas[i], x) * ang_bracket(lambdas[i], y)
        if c_val == 0:
            raise ValueError(f"Reference spinor orthogonal to particle {i}")
        C.append(c_val)
        
    # 3. Build Weighted Laplacian
    L_tilde = matrix(QQ, n, n)
    
    def get_w(i, j):
        if i > j: return weights.get((j, i), 0)
        return weights.get((i, j), 0)
        
    for i in range(n):
        row_sum = 0
        for j in range(n):
            if i != j:
                w = get_w(i, j)
                val = w * C[i] * C[j]
                L_tilde[i, j] = -val
                row_sum += val
        L_tilde[i, i] = row_sum
        
    # 4. Build Phi_tilde = -D^-1 L_tilde D^-1
    # (Phi_tilde)_ij = - (1/Ci) (L_tilde)_ij (1/Cj)
    #                = - (1/Ci) (-w_ij Ci Cj) (1/Cj)
    #                = w_ij  (for i != j)
    # (Phi_tilde)_ii = - (1/Ci^2) L_tilde_ii
    #                = - (1/Ci^2) sum (w_ik Ci Ck)
    #                = - sum (w_ik Ck/Ci)
    
    Phi_tilde = matrix(QQ, n, n)
    for i in range(n):
        for j in range(n):
            Phi_tilde[i, j] = - L_tilde[i, j] / (C[i] * C[j])
            
    return L_tilde, C, Phi_tilde


def hodges_minor_matrix_tree(lambdas, tilde_lambdas):
    """
    Legacy wrapper for 6-point tree sum using enumeration.
    Kept for compatibility with Phase D scripts.
    """
    n = 6
    weights = {}
    for i in range(n):
        for j in range(i + 1, n):
            ang = ang_bracket(lambdas[i], lambdas[j])
            sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
            weights[(i, j)] = sq / ang
            
    return tree_sum_kirchhoff(weights, n, delete=0)

if __name__ == "__main__":
    # Quick sanity check
    print("Running matrix_tree.py sanity check...")
    pass
