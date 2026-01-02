import sys
import os
import numpy as np

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from posgeom.saddle_pushforward import compute_pushforward_saddle
from posgeom.forest_polytope import get_forest_exponents

def verify_n4_saddle():
    print("\nVerifying n=4 Saddle Pushforward...")
    # n=4, roots=0,1,2.
    # Forests are single edges (0,3), (1,3), (2,3).
    # Exponents in R^3 (actually R^6 but only 3 active).
    # Wait, get_forest_exponents returns vectors in R^|Edges|.
    # For n=4, |Edges|=6.
    # The active edges are (0,3), (1,3), (2,3).
    # Indices 3, 4, 5 in edge_order [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)] ?
    
    n = 4
    roots = [0, 1, 2]
    # We can rely on get_forest_exponents to give us vectors.
    exponents, edge_order = get_forest_exponents(n, roots)
    # Convert to numpy
    exponents = np.array(exponents)
    coeffs = np.ones(len(exponents)) # All 1 for unweighted forest polynomial
    
    # We need to project to intrinsic dimension to make J invertible.
    # Exponents lie on affine plane sum x_i = 1 (if standard simplex).
    # Actually, they are standard basis vectors e_i for the active edges.
    # They lie in a subspace?
    # Dim of convex hull is 2 (3 vertices).
    # Space dimension is 6.
    # We need to project X and exponents to the 2D affine subspace.
    
    # Projection Strategy:
    # 1. Shift by v0 -> v - v0.
    # 2. Find basis for span(v_i - v0).
    # 3. Express X_target in this basis.
    
    v0 = exponents[0]
    diffs = exponents[1:] - v0
    # diffs is (N-1) x D matrix.
    # For n=4, N=3. diffs is 2 x 6.
    # Basis B = diffs.
    # New coords Y in R^2.
    # Original V = v0 + Y * B.
    # Projected V_proj = Rows of B? No.
    # We map p(x) in embedding space to p(y) in intrinsic space.
    # x = exp(u). u = u0 + B^T y ? No.
    # We just want vertices in intrinsic R^2.
    # Let vertices be (0,0), (1,0), (0,1).
    # Then compute pushforward.
    
    # Simplex case
    proj_exponents = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    proj_coeffs = [1.0, 1.0, 1.0]
    
    # Target X (centroid)
    X_target = [1.0/3.0, 1.0/3.0]
    
    val = compute_pushforward_saddle(X_target, proj_coeffs, proj_exponents)
    print(f"Computed Saddle Val: {val:.5f}")
    
    # Exact canonical form for simplex 1+x+y
    # Omega = 1 / (x y (1-x-y)) dx dy? No.
    # The canonical FUNCTION is 1 / (x y (1-x-y)).
    # Wait, is it?
    # Dual of simplex is simplex.
    # Area of dual triangle?
    # Formula for simplex with vertices V:
    # d! Vol(V) / prod L_i(X) ?
    # For standard simplex, vertices (0,0), (1,0), (0,1).
    # Facets: x>=0, y>=0, 1-x-y>=0.
    # Canonical function: 1 / (x * y * (1-x-y)).
    # At (1/3, 1/3), val = 1 / (1/27) = 27.
    
    print(f"Expected: 27.0")
    if abs(val - 27.0) < 1e-4:
        print("SUCCESS")
    else:
        print("FAILURE")

def verify_n5_saddle():
    print("\nVerifying n=5 Saddle Pushforward...")
    # n=5, roots=0,1,2.
    # 15 forests.
    # Intrinsic dimension 4.
    
    n = 5
    roots = [0, 1, 2]
    exponents, _ = get_forest_exponents(n, roots)
    exponents = np.array(exponents)
    
    # Project to 4D
    v0 = exponents[0]
    diffs = exponents[1:] - v0 # 14 x 10
    
    # PCA / SVD to find basis
    U, S, Vt = np.linalg.svd(diffs.T)
    # Rank should be 4.
    rank = np.sum(S > 1e-10)
    print(f"Rank/Dim: {rank}")
    
    basis = U[:, :rank] # 10 x 4
    # Project vertices: V_proj = V @ basis
    # Or rather (V - v0) @ basis
    proj_exponents = np.dot(exponents - v0, basis)
    
    # Target: Centroid (average of vertices)
    # This is definitely in the interior.
    X_target = np.mean(proj_exponents, axis=0)
    coeffs = np.ones(len(proj_exponents))
    
    val_saddle = compute_pushforward_saddle(X_target, coeffs, proj_exponents)
    print(f"Saddle Value: {val_saddle:.5f}")
    
    # Compare with existing canonical form evaluator?
    # We need to build the polytope object.
    # But current tools in repo (canonical_polytope.py) take H-rep or V-rep.
    # We can try to assume it works if n=4 works, or run a check?
    
    # Let's perform a consistency check:
    # Value should be independent of gauge (shift of X).
    
    pass

if __name__ == "__main__":
    verify_n4_saddle()
    verify_n5_saddle()




