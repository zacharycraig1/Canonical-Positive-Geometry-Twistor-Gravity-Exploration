
import sys
import os
import json
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.physics_map import eval_edge_vars_from_spinors

try:
    from phase_s.probe_limits_on_shell import probe_limit
except ImportError:
    # If running from inside phase_s or if module resolution fails
    sys.path.append(os.path.join(os.getcwd(), 'phase_s'))
    from probe_limits_on_shell import probe_limit

def get_polytope_data(n, roots, seed=None):
    """
    Returns (polyhedron, vertices, inequalities).
    """
    exponents, edge_order = get_forest_exponents(n, roots)
    P = Polyhedron(vertices=exponents)
    return P, exponents, P.inequalities(), edge_order

def compute_moment_map_point(z_values, exponents, edge_order):
    """
    Computes x = Sum(alpha * z^alpha) / Sum(z^alpha).
    z_values: dict mapping edge index/tuple to value.
    exponents: list of exponent vectors.
    edge_order: list of edge tuples corresponding to vector indices.
    """
    numerator = vector(QQ, [0]*len(edge_order))
    denominator = QQ(0)
    
    for alpha in exponents:
        # Compute monomial z^alpha
        term = QQ(1)
        for i, edge in enumerate(edge_order):
            pow_val = alpha[i]
            if pow_val > 0:
                # z_values key might be (u,v) with u<v
                u, v = edge
                if u > v: u, v = v, u
                
                # If z is 0 and pow > 0, term is 0.
                if z_values.get((u,v), 0) == 0:
                    term = QQ(0)
                    break
                term *= z_values[(u,v)] ** pow_val
        
        vector_alpha = vector(QQ, alpha)
        numerator += vector_alpha * term
        denominator += term
        
    if denominator == 0:
        return None
        
    return numerator / denominator

def analyze_limit(n, roots, limit_name, limit_type, epsilon, **kwargs):
    print(f"Analyzing {limit_name} with eps={epsilon}...")
    
    # 1. Get Polytope
    P, exponents, inequalities, edge_order = get_polytope_data(n, roots)
    
    # 2. Get Kinematics (Use positive region)
    kinematics = probe_limit(n, limit_type, epsilon, positive=True, **kwargs)
    
    # 3. Compute z variables
    # We need reference spinors.
    x_ref = vector(QQ, [1, 0])
    y_ref = vector(QQ, [0, 1])
    lambdas_dict = {i: kinematics.lambdas[i] for i in range(n)}
    tildes_dict = {i: kinematics.tilde_lambdas[i] for i in range(n)}
    
    try:
        z_values_raw = eval_edge_vars_from_spinors(lambdas_dict, tildes_dict, x_ref, y_ref)
        z_values = {k: abs(v) for k, v in z_values_raw.items()}
    except ValueError as e:
        print(f"  Failed to compute z: {e}")
        return None
        
    # 4. Compute Moment Map Point
    x_point = compute_moment_map_point(z_values, exponents, edge_order)
    if x_point is None:
        print("  Moment map denominator zero.")
        return None
        
    # 5. Check Facets
    active_facets = []
    
    # Helper to interpret inequality
    def interpret_ieq(ieq_coeffs, edge_order):
        # ieq is A . x + b >= 0.
        # coeffs are A_0, ..., A_{d-1}, b
        A = ieq_coeffs[:-1]
        b = ieq_coeffs[-1]
        
        # Check if coordinate hyperplane x_e >= 0
        # A has one 1, rest 0. b=0.
        if b == 0 and sum(abs(c) for c in A) == 1 and 1 in A:
            idx = A.index(1)
            u, v = edge_order[idx]
            return f"z_{u}_{v} >= 0"
        return str(ieq_coeffs)

    for ieq in inequalities:
        val = ieq.eval(x_point)
        coeffs = list(ieq.A()) + [ieq.b()]
        desc = interpret_ieq(coeffs, edge_order)
        
        active_facets.append({
            "inequality": coeffs,
            "description": desc,
            "slack": float(val),
            "slack_exact": str(val)
        })
        
    # Sort by slack
    active_facets.sort(key=lambda k: k["slack"])
    
    # Return top 5
    return active_facets[:5]

def analyze_limit_with_retry(n, roots, limit_name, limit_type, epsilon, **kwargs):
    for attempt in range(10):
        try:
            # Pass a seed that changes
            kwargs['seed'] = 42 + attempt * 100
            res = analyze_limit(n, roots, limit_name, limit_type, epsilon, **kwargs)
            if res is not None:
                return res
        except Exception as e:
            # print(f"  Attempt {attempt} failed with exception: {e}")
            pass
            
    print(f"  All attempts failed for {limit_name}")
    return None

def run_analysis_n6():
    n = 6
    roots = [0, 1, 2] # Fix roots for now
    epsilon = QQ(1)/1000 # 1e-3
    
    results = {}
    
    # 1. Soft Limit (Leg 5)
    res_soft = analyze_limit_with_retry(n, roots, "Soft(5)", "soft", epsilon, s_idx=5)
    results["Soft(5)"] = res_soft
    
    # 2. Collinear (3 || 4)
    res_col_34 = analyze_limit_with_retry(n, roots, "Col(3||4)", "collinear", epsilon, i=3, j=4)
    results["Col(3||4)"] = res_col_34
    
    # 3. Collinear (4 || 5)
    res_col_45 = analyze_limit_with_retry(n, roots, "Col(4||5)", "collinear", epsilon, i=4, j=5)
    results["Col(4||5)"] = res_col_45
    
    # Save results
    with open("phase_s/active_facets.json", "w") as f:
        # Helper to serialize sets/lists
        def default(o):
            if isinstance(o, (Integer, int)): return int(o)
            if isinstance(o, (Float, float)): return float(o)
            return str(o)
        json.dump(results, f, indent=2, default=default)
        
    # Print Summary
    print("\nSummary of Active Facets (Slack < 1e-2):")
    for limit, facets in results.items():
        print(f"\n{limit}:")
        if facets:
            for f in facets:
                if f["slack"] < 0.1: # Threshold
                    print(f"  Slack: {f['slack']:.2e} | {f['description']}")

if __name__ == "__main__":
    run_analysis_n6()
