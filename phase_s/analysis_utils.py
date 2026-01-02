
import sys
import os
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.physics_map import eval_edge_vars_from_spinors

try:
    from phase_s.probe_limits_on_shell import probe_limit
except ImportError:
    sys.path.append(os.path.join(os.getcwd(), 'phase_s'))
    from probe_limits_on_shell import probe_limit

def get_polytope_data(n, roots):
    """
    Returns (polyhedron, vertices, inequalities, edge_order).
    """
    exponents, edge_order = get_forest_exponents(n, roots)
    P = Polyhedron(vertices=exponents)
    return P, exponents, P.inequalities(), edge_order

def compute_moment_map_point(z_values, exponents, edge_order):
    """
    Computes x = Sum(alpha * z^alpha) / Sum(z^alpha).
    """
    numerator = vector(QQ, [0]*len(edge_order))
    denominator = QQ(0)
    
    for alpha in exponents:
        # Compute monomial z^alpha
        term = QQ(1)
        for i, edge in enumerate(edge_order):
            pow_val = alpha[i]
            if pow_val > 0:
                u, v = edge
                if u > v: u, v = v, u
                
                val = z_values.get((u,v), 0)
                if val == 0:
                    term = QQ(0)
                    break
                term *= val ** pow_val
        
        vector_alpha = vector(QQ, alpha)
        numerator += vector_alpha * term
        denominator += term
        
    if denominator == 0:
        return None
        
    return numerator / denominator

def interpret_ieq(ieq_coeffs, edge_order):
    A = ieq_coeffs[:-1]
    b = ieq_coeffs[-1]
    
    # Check if coordinate hyperplane x_e >= 0
    if b == 0 and sum(abs(c) for c in A) == 1 and 1 in A:
        idx = A.index(1)
        u, v = edge_order[idx]
        return f"z_{u}_{v} >= 0"
    return str(ieq_coeffs)

def analyze_limit(n, roots, limit_name, limit_type, epsilon, **kwargs):
    # 1. Get Polytope
    P, exponents, inequalities, edge_order = get_polytope_data(n, roots)
    
    # 2. Get Kinematics (Use positive region)
    try:
        kinematics = probe_limit(n, limit_type, epsilon, positive=True, **kwargs)
    except Exception as e:
        # print(f"  Probe failed: {e}")
        return None
    
    # 3. Compute z variables
    x_ref = vector(QQ, [1, 0])
    y_ref = vector(QQ, [0, 1])
    lambdas_dict = {i: kinematics.lambdas[i] for i in range(n)}
    tildes_dict = {i: kinematics.tilde_lambdas[i] for i in range(n)}
    
    try:
        z_values_raw = eval_edge_vars_from_spinors(lambdas_dict, tildes_dict, x_ref, y_ref)
        z_values = {k: abs(v) for k, v in z_values_raw.items()}
    except ValueError as e:
        return None
        
    # 4. Compute Moment Map Point
    x_point = compute_moment_map_point(z_values, exponents, edge_order)
    if x_point is None:
        return None
        
    # 5. Check Facets
    active_facets = []
    
    for ieq in inequalities:
        val = ieq.eval(x_point)
        coeffs = list(ieq.A()) + [ieq.b()]
        desc = interpret_ieq(coeffs, edge_order)
        
        active_facets.append({
            "inequality": coeffs,
            "description": desc,
            "slack": float(val)
        })
        
    # Sort by slack
    active_facets.sort(key=lambda k: k["slack"])
    
    # Return top 5
    return active_facets[:5]

def analyze_limit_with_retry(n, roots, limit_name, limit_type, epsilon, **kwargs):
    for attempt in range(5):
        try:
            kwargs['seed'] = 42 + attempt * 100
            res = analyze_limit(n, roots, limit_name, limit_type, epsilon, **kwargs)
            if res is not None:
                return res
        except Exception:
            pass
            
    return None




