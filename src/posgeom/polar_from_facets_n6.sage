import sys
import os
import json
from sage.all import QQ, matrix, vector, MixedIntegerLinearProgram, Polyhedron

sys.path.append(os.getcwd())

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)

def construct_polar():
    print("Loading hull data...")
    hull_eq = load_json("RESULTS/facets_n6_eq_exact.json")
    facets = load_json("RESULTS/facet_dictionary_n6.json")
    
    # We are in 15 dim, but polytope is 11 dim.
    # To construct polar, we should project to 11 dim intrinsic space?
    # Or work in 15 dim but polar will be infinite in normal directions?
    # Standard polar duality is for full-dim polytopes.
    # We should project to the affine hull.
    
    # Parametrization from hull equations
    # X = X0 + N * t
    # We found X0 and N in fit script. Let's re-derive or load.
    
    equations = hull_eq["equations"]
    # ... Solve E X = -b ...
    # Easier to just load the fit result which has X0, N.
    if os.path.exists("RESULTS/kinematic_map_n6_fit.json"):
        fit = load_json("RESULTS/kinematic_map_n6_fit.json")
        X0 = vector(QQ, fit["X0"])
        # N_rows contains the 15 rows of the 15x11 matrix N
        N = matrix(QQ, fit["N_rows"]) # 15x11
    else:
        print("Fit map not found, cannot define subspace.")
        return

    print("Projecting facets to intrinsic t-space (dim 11)...")
    # L_F(X) = b + A X = b + A(X0 + N t) = (b + A X0) + (A N) t
    # Let beta_F = b + A X0
    # Let alpha_F = A N (row vector)
    # Inequality: beta_F + alpha_F . t >= 0
    
    intrinsic_inequalities = [] # Format for Polyhedron: [beta, alpha...] so beta + alpha.x >= 0
    
    for f in facets:
        f_coeffs = [QQ(x) for x in f["ineq_b_A"]]
        b_F = f_coeffs[0]
        A_F = vector(QQ, f_coeffs[1:])
        
        beta = b_F + A_F.dot_product(X0)
        alpha = A_F * N # vector length 11
        
        # Sage Polyhedron expects [b, a1, ... an] for b + A x >= 0
        ieq = [beta] + list(alpha)
        intrinsic_inequalities.append(ieq)
        
    print(f"Constructing Polyhedron from {len(intrinsic_inequalities)} inequalities...")
    P = Polyhedron(ieqs=intrinsic_inequalities, base_ring=QQ)
    
    if P.is_empty():
        print("Polyhedron is empty!")
        return
        
    print(f"Polyhedron constructed. Vertices: {P.n_vertices()}")
    
    # Find interior point (Chebyshev center)
    # Sage doesn't have direct chebyshev center?
    # We can use P.center() if bounded?
    # Or just average of vertices.
    
    # Vertices might be many?
    if P.n_vertices() > 10000:
        print("Too many vertices to compute center from all.")
        center = P.vertices()[0].vector() # unsafe
    else:
        verts = P.vertices()
        center = sum(v.vector() for v in verts) / len(verts)
        
    print(f"Center t_c: {center}")
    
    # Check if center is strictly interior
    min_slack = min(beta + vector(QQ, alpha).dot_product(center) for beta, *alpha in intrinsic_inequalities)
    print(f"Min slack at center: {min_slack}")
    
    if min_slack <= 0:
        print("Center is not strictly interior. Using LP to find interior point.")
        # Setup LP
        p = MixedIntegerLinearProgram(maximization=True)
        t = p.new_variable(real=True)
        r = p.new_variable(real=True, nonnegative=True)
        p.set_objective(r[0])
        
        # For each ineq: beta + alpha.t >= r
        for ieq in intrinsic_inequalities:
            beta = ieq[0]
            alpha = ieq[1:]
            p.add_constraint(beta + sum(alpha[i]*t[i] for i in range(11)) >= r[0])
            
        try:
            p.solve()
            center = vector(QQ, p.get_values(t))
            min_slack = p.get_values(r)[0]
            print(f"LP Center found. Slack: {min_slack}")
        except Exception as e:
            print(f"LP failed: {e}")
            return

    # Shift to center
    # t = center + y
    # Ineq: beta + alpha(center + y) >= 0
    # (beta + alpha.center) + alpha.y >= 0
    # Let const_new = beta + alpha.center (positive)
    # Normalize: 1 + (alpha / const_new) . y >= 0
    # 1 >= - (alpha / const_new) . y
    # 1 >= w . y  where w = - alpha / const_new
    # This is standard form a.y <= 1.
    
    polar_vertices = []
    
    for ieq in intrinsic_inequalities:
        beta = ieq[0]
        alpha = vector(QQ, ieq[1:])
        const_val = beta + alpha.dot_product(center)
        
        if abs(const_val) < 1e-10:
            # Facet passes through center?? Should not happen if slack > 0
            continue
            
        w = - alpha / const_val
        polar_vertices.append(list(w))
        
    print(f"Computed {len(polar_vertices)} polar vertices.")
    
    # Save polar vertices
    # Also save center and N map to reconstruct full X
    # X_center = X0 + N * center
    
    def json_friendly(obj):
        if hasattr(obj, "numerator"):
            return float(obj)
        if isinstance(obj, list):
            return [json_friendly(x) for x in obj]
        return obj

    result = {
        "polar_vertices": json_friendly(polar_vertices),
        "center_t": json_friendly(list(center)),
        "center_slack": float(min_slack)
    }
    
    with open("RESULTS/polar_polytope_n6.json", "w") as f:
        json.dump(result, f, indent=2)
    print("Saved RESULTS/polar_polytope_n6.json")

if __name__ == "__main__":
    construct_polar()

