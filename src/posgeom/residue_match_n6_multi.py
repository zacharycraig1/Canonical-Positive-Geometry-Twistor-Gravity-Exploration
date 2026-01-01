import sys
import os
import json
import argparse
import random
import math
from sage.all import RR, vector, matrix, QQ

sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.intrinsic_lattice import IntrinsicLattice
from src.posgeom.moment_map_laplacian import MomentMapLaplacian
from src.posgeom.physics_map import eval_edge_vars_from_spinors

def load_facets():
    path = "RESULTS/facet_dictionary_n6.json"
    if not os.path.exists(path):
        print(f"Warning: {path} not found. Using empty facet list.")
        return []
    with open(path, "r") as f:
        return json.load(f)

def compute_slacks(X_vector, facets_data):
    """
    Computes slack u_F = b + A.X for each facet.
    facets_data is list of dicts with "ineq_b_A".
    """
    slacks = {}
    for entry in facets_data:
        fid = entry["facet_id"]
        coeffs = [float(x) for x in entry["ineq_b_A"]] # [b, A0, A1...]
        b = coeffs[0]
        A = vector(RR, coeffs[1:])
        
        # X is vector(QQ) or RR
        # Ensure dimensions match
        if len(A) != len(X_vector):
            # Might happen if X is ambient vs intrinsic?
            # The inequalities are usually in AMBIENT space z variables.
            # X from MomentMap is in ambient space (log z_e basis).
            pass
            
        slack = b + A.dot_product(vector(RR, X_vector))
        slacks[fid] = slack
    return slacks

def run_residue_check_multi():
    n = 6
    roots = [0, 1, 2]
    
    print("Loading geometry...")
    exponents, edge_order = get_forest_exponents(n, roots)
    lattice = IntrinsicLattice(exponents)
    mml = MomentMapLaplacian(n, roots, edge_order)
    facets_data = load_facets()
    
    # Probes: (name, epsilon_func)
    probes = []
    
    # 1. Collinear (0,1) -> s_01 -> 0
    def probe_01(eps):
        ts = [1, 2, 3, 4, 5, 6] # Base
        ts[1] = ts[0] + eps
        return ts
    probes.append(("Collinear 0,1", probe_01))
    
    # 2. Collinear (3,4) (Non-roots) -> s_34 -> 0
    def probe_34(eps):
        ts = [1, 2, 3, 4, 5, 6]
        ts[4] = ts[3] + eps
        return ts
    probes.append(("Collinear 3,4", probe_34))
    
    # 3. Triple (0,1,2) -> s_012 -> 0
    def probe_012(eps):
        ts = [1, 2, 3, 4, 5, 6]
        ts[1] = ts[0] + eps
        ts[2] = ts[0] + 2*eps
        return ts
    probes.append(("Triple 0,1,2", probe_012))
    
    # 4. Triple (3,4,5) -> s_345 -> 0
    def probe_345(eps):
        ts = [1, 2, 3, 4, 5, 6]
        ts[4] = ts[3] + eps
        ts[5] = ts[3] + 2*eps
        return ts
    probes.append(("Triple 3,4,5", probe_345))

    epsilons = [1e-2, 1e-3, 1e-4, 1e-5]
    
    # Reference spinors
    x = vector(RR, [1, -2])
    y = vector(RR, [1, 12])
    
    for name, func in probes:
        print(f"\n--- Probe: {name} ---")
        
        results = []
        
        for eps in epsilons:
            ts = func(eps)
            lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
            # Tildes random but fixed per eps? Or just generic?
            # For 1/s^2 scaling in gravity, we need momentum conservation.
            # But here we are using Twistor variables, so momenta are derived.
            # M_MHV is computed from Hodges which is valid for these lambdas/tildes.
            # Just need tildes to be generic.
            random.seed(42) # Consistent tildes
            tildes = {i: vector(RR, [1, random.uniform(0, 10)]) for i in range(n)}
            
            # Compute z
            try:
                z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
                z_vals = [z_map[(u,v)] for u,v in edge_order]
            except ValueError:
                continue
                
            # Compute X, H, Omega
            try:
                X, H = mml.compute_X_H(z_vals)
                B_mat = lattice.B
                H_int = B_mat.transpose() * H * B_mat
                det_H = H_int.det()
                Omega = 1.0/det_H if abs(det_H) > 1e-20 else float('inf')
                
                # Slacks
                slacks = compute_slacks(X, facets_data)
                min_slack = min(slacks.values()) if slacks else 0
                
                # Find facets approaching 0 (slack < 1.0 roughly? No, < 0.1)
                near_facets = [fid for fid, s in slacks.items() if abs(s) < 0.1]
                
                results.append({
                    "eps": eps,
                    "Omega": Omega,
                    "min_slack": min_slack,
                    "near_facets": near_facets
                })
                
                print(f"eps={eps:.1e} | Omega={Omega:.2e} | MinSlack={min_slack:.2e} | NearFacets={near_facets}")
                
            except Exception as e:
                print(f"Error at eps={eps}: {e}")
                
        # Analyze scaling
        if len(results) >= 2:
            r1 = results[0]
            r2 = results[-1]
            if r1["Omega"] != 0 and r2["Omega"] != 0:
                slope_Om = (math.log(abs(r2["Omega"])) - math.log(abs(r1["Omega"]))) / (math.log(r2["eps"]) - math.log(r1["eps"]))
            else:
                slope_Om = 0
                
            if r1["min_slack"] != 0 and r2["min_slack"] != 0:
                slope_Slack = (math.log(abs(r2["min_slack"])) - math.log(abs(r1["min_slack"]))) / (math.log(r2["eps"]) - math.log(r1["eps"]))
            else:
                slope_Slack = 0
                
            print(f"Scaling: Omega ~ eps^{slope_Om:.2f}, MinSlack ~ eps^{slope_Slack:.2f}")
            
            if slope_Slack > 0.5:
                print("SUCCESS: Facet approached!")
            else:
                print("FAILURE: No facet approach detected (MinSlack constant).")

if __name__ == "__main__":
    run_residue_check_multi()



