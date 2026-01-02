
import sys
import os
import json
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.canonical_form import CanonicalFormEvaluator
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from phase_s.probe_limits_on_shell import probe_limit

def check_residue_col45():
    print("S2: Checking Iterated Residue for Col(4||5)...")
    
    # 1. Load active facets
    with open("phase_s/active_facets.json", "r") as f:
        data = json.load(f)
        
    facets_data = data.get("Col(4||5)", [])
    if not facets_data:
        print("No active facets found for Col(4||5)")
        return
        
    print(f"Found {len(facets_data)} active inequalities.")
    
    # 2. Build Polyhedron (n=6, roots=[0,1,2])
    n = 6
    roots = [0, 1, 2]
    exponents, edge_order = get_forest_exponents(n, roots)
    P = Polyhedron(vertices=exponents)
    
    # 3. Identify Face Vertices
    # A vertex v is on the face if it satisfies A.v + b = 0 for all active facets
    # Note: JSON stores slack. We want slack=0 (or epsilon).
    
    active_ineqs = []
    for entry in facets_data:
        if entry['slack'] < 1e-5: # Threshold
            active_ineqs.append(entry['inequality'])
            
    print(f"Number of binding inequalities: {len(active_ineqs)}")
    
    face_vertices = []
    for v in exponents:
        v_vec = vector(QQ, v)
        on_face = True
        for coeffs in active_ineqs:
            # coeffs is [A0...Ad, b]
            # ineq is A.x + b >= 0
            A = vector(QQ, coeffs[:-1])
            b = QQ(coeffs[-1])
            val = A.dot_product(v_vec) + b
            if abs(val) > 1e-10:
                on_face = False
                break
        if on_face:
            face_vertices.append(v)
            
    print(f"Face has {len(face_vertices)} vertices out of {len(exponents)}.")
    
    if len(face_vertices) == 0:
        print("Error: Empty face!")
        return
        
    F = Polyhedron(vertices=face_vertices)
    print(f"Face Dimension: {F.dim()}")
    print(f"Polytope Dimension: {P.dim()}")
    print(f"Codimension: {P.dim() - F.dim()}")
    
    # 4. Evaluate Omega(F) at a point on the face
    # We need a W that lies in the face?
    # No, W is dual.
    # If we are verifying the residue, we evaluate Omega(F) at the limit point projected to F?
    # The physics map sends the collinear limit to this face.
    # So we take the limit point z(epsilon -> 0).
    
    # Get limit point
    limit_res = probe_limit(n, "collinear", epsilon=1e-8, i=4, j=5, positive=True, seed=42)
    lambdas = limit_res.lambdas
    tildes = limit_res.tilde_lambdas
    
    # Compute z
    x_ref = vector(QQ, [1, 0])
    y_ref = vector(QQ, [0, 1])
    z_vals = eval_edge_vars_from_spinors(lambdas, tildes, x_ref, y_ref)
    
    # Construct W vector from z_vals
    W_comps = []
    # W must match edge order
    # Note: W is homogeneous dual.
    # z variables correspond to coordinate planes?
    # Actually, the map is z_e.
    # The point in the dual space is Y = (z_e).
    # Wait, my CanonicalFormEvaluator expects W.
    # If the polytope is in the space of exponents, W are the log-variables z?
    # No. The polytope P is the Newton polytope.
    # The form is Omega(Y) where Y are the variables z_e.
    # Yes. So W = vector(z_e).
    # But W must be homogeneous?
    # If P is in R^E, then W is in (P^E)*?
    # Actually, usually W = (1, z_e).
    # Let's assume W = (1, z_e) for now.
    
    W_list = [1] # Homogeneous component
    for u, v in edge_order:
        val = z_vals.get((u, v))
        if val is None: val = z_vals.get(f"z_{u}_{v}")
        W_list.append(val)
        
    W = vector(QQ, W_list)
    
    # Evaluate Omega(F) at W
    # Note: Omega(F) is defined on the subspace spanned by F.
    # The Evaluator handles projection!
    
    print("Evaluating Omega(F) at limit point...")
    try:
        val = CanonicalFormEvaluator.eval_polytope(F, W)
        print(f"Omega(F) value: {val}")
        if val != 0:
            print("  [PASS] Non-zero residue.")
        else:
            print("  [WARN] Residue is zero?")
    except Exception as e:
        print(f"  [FAIL] Evaluation error: {e}")

if __name__ == "__main__":
    check_residue_col45()




