import json
import os
import sys
import numpy as np
from sage.all import *

sys.path.append(os.getcwd())
# Import physics map or reimplement
# z_ij = [ij]/<ij> * C_i * C_j
# We need to define C_i carefully. Reference spinors?
# src/posgeom/physics_map.py likely has it.
from src.posgeom.physics_map import eval_edge_vars_from_spinors

def generate_random_spinors(n=6):
    # Random complex spinors
    # Convert to dict for compatibility if needed, but eval_edge_vars uses indexing [i]
    lam = [vector(CC, [random(), random()]) for _ in range(n)]
    tlam = [vector(CC, [random(), random()]) for _ in range(n)]
    return lam, tlam

def apply_soft_limit(lam, tlam, leg_idx, epsilon):
    # Scale leg p -> epsilon^2 p
    # lam -> epsilon * lam
    # tlam -> epsilon * tlam
    n = len(lam)
    new_lam = [v if i != leg_idx else v * epsilon for i, v in enumerate(lam)]
    new_tlam = [v if i != leg_idx else v * epsilon for i, v in enumerate(tlam)]
    return new_lam, new_tlam

def apply_collinear_limit(lam, tlam, i, j, epsilon):
    # lam_i || lam_j
    # Set lam_i = lam_j + epsilon * error
    # tlam must adjust to conserve momentum? 
    # For "probe", we don't strictly need momentum conservation if we just map spinors -> z.
    # The map z(spinors) is defined for any spinors (off-shell).
    # Momentum conservation imposes constraints ON z.
    # But the pushforward geometry should capture the limit BEHAVIOR.
    
    # Simple collinear:
    # lam_i = lam_j
    # But we want approach.
    new_lam = list(lam)
    # Make i close to j
    new_lam[i] = new_lam[j] + vector(CC, [random(), random()]) * epsilon
    
    # Keep tlam generic?
    return new_lam, tlam

def apply_bcfw_limit(lam, tlam, i, j, epsilon):
    # BCFW shift sends a propagator to 0.
    # Usually s_S -> 0.
    # We want to probe a factorization channel.
    # e.g. s_{123} -> 0.
    # This requires specific solving.
    # Easier: Just manually set P_{123}^2 = epsilon.
    # Hard to do with just spinors directly without solving.
    
    # Alternative: Use "subset limit" where we scale a subset of particles?
    # No, that's multi-soft.
    
    # Let's stick to Soft and Collinear for now as they are primary boundaries.
    return None, None

def analyze_scaling(val1, val2, eps1, eps2):
    # val ~ C * eps^alpha
    # log(val) ~ log(C) + alpha * log(eps)
    # alpha = (log(val1) - log(val2)) / (log(eps1) - log(eps2))
    
    if val1 == 0 or val2 == 0:
        return 0 # Or -inf?
    
    # Use magnitude
    lv1 = float(log(abs(val1)))
    lv2 = float(log(abs(val2)))
    le1 = float(log(eps1))
    le2 = float(log(eps2))
    
    alpha = (lv1 - lv2) / (le1 - le2)
    return alpha

def run_probes():
    n = 6
    roots = [0, 1, 2] # Convention
    
    # Load boundary dict to match results
    dict_path = os.path.join("RESULTS", f"boundary_dictionary_n{n}.json")
    with open(dict_path, 'r') as f:
        boundary_dict = json.load(f)
        
    # Debug: list available upper bounds
    print("Available Upper Bound Facets:")
    for fid, info in boundary_dict.items():
        if info['type'] == 'upper_bound':
            print(f"  Facet {fid}: {info['edge']}")

    probes = []
    
    # 1. Soft Limit (Leg 5)
    print("Probing Soft Limit (Leg 5)...")
    lam0, tlam0 = generate_random_spinors(n)
    
    eps1, eps2 = 1e-5, 1e-8
    
    # Soft leg 5
    l1, tl1 = apply_soft_limit(lam0, tlam0, 5, eps1)
    l2, tl2 = apply_soft_limit(lam0, tlam0, 5, eps2)
    
    # Compute z
    # We need auxiliary spinors x, y for C_i. 
    # physics_map usually handles this or we define them.
    # Let's inspect src/posgeom/physics_map.py to see API.
    # Assuming compute_z_variables_from_spinors(lam, tlam, x=None, y=None)
    
    # Mock reference spinors
    x = vector(CC, [1, 0])
    y = vector(CC, [0, 1])
    
    z1 = eval_edge_vars_from_spinors(l1, tl1, x, y)
    z2 = eval_edge_vars_from_spinors(l2, tl2, x, y)
    
    # Helper to map z dict to edge indices
    def get_exponents(z_dict_1, z_dict_2):
        exps = {}
        for k in z_dict_1:
            if isinstance(k, str): continue # Skip string keys
            alpha = analyze_scaling(z_dict_1[k], z_dict_2[k], eps1, eps2)
            exps[k] = alpha
        return exps
        
    exps_soft = get_exponents(z1, z2)
    
    probes.append({
        "name": "Soft Leg 5",
        "exponents": exps_soft
    })
    
    # 2. Collinear Limit (3 || 4)
    print("Probing Collinear Limit (3 || 4)...")
    l1, tl1 = apply_collinear_limit(lam0, tlam0, 3, 4, eps1)
    l2, tl2 = apply_collinear_limit(lam0, tlam0, 3, 4, eps2)
    
    z1 = eval_edge_vars_from_spinors(l1, tl1, x, y)
    z2 = eval_edge_vars_from_spinors(l2, tl2, x, y)
    
    exps_coll = get_exponents(z1, z2)
    
    probes.append({
        "name": "Collinear 3||4",
        "exponents": exps_coll
    })
    
    # Match to facets
    # For each probe, find facets that are "triggered".
    # A facet is triggered if the variables involved in the inequality scale significantly?
    # Or if the forest polynomial on that facet dominates?
    
    # Better: Check if z satisfies the facet inequality "approaching boundary".
    # e.g. lower bound x_ij >= 0. If z_ij -> 0, then we are at x_ij = 0?
    # Map: z_ij = e^{x_ij}? Or similar.
    # If z -> 0, x -> -inf ??
    # The map from Polytope to Kinematics is \Phi.
    # \Phi_*(\Omega_P) = \Omega_{phys}.
    # Usually this means z_ij corresponds to edge weights.
    # Small z_ij corresponds to x_ij boundaries?
    # If z_ij are Arkani-Hamed/Trnka variables:
    # z_ij small -> x_ij boundary?
    
    # Let's just record which z_ij go to 0 or infinity.
    # And see which facets involve those edges.
    
    results = []
    
    for probe in probes:
        print(f"Analyzing {probe['name']}...")
        scaling_vars = {k: v for k, v in probe['exponents'].items() if abs(v) > 0.1}
        
        # Identify "Vanishing" (alpha > 0) and "Diverging" (alpha < 0)
        vanishing = [k for k, v in scaling_vars.items() if v > 0.5]
        diverging = [k for k, v in scaling_vars.items() if v < -0.5]
        
        print(f"  Vanishing: {vanishing}")
        print(f"  Diverging: {diverging}")
        
        # Check against dictionary
        matched_facets = []
        for fid, info in boundary_dict.items():
            # Check if this facet's definition matches the behavior
            is_match = False
            
            if info['type'] == 'lower_bound':
                # Expect corresponding edge to VANISH (z -> 0)
                e = tuple(sorted(info['edge']))
                if e in vanishing:
                    is_match = True
                    
            elif info['type'] == 'upper_bound':
                # Expect corresponding edge to DIVERGE (z -> inf)
                e = tuple(sorted(info['edge']))
                if e in diverging:
                    is_match = True
            
            elif info['type'] == 'subset_rank':
                # Subset S.
                # If sum x_e <= rank -> z_e small?
                # Usually corresponds to multi-particle channel.
                # We haven't probed BCFW yet.
                pass
                
            if is_match:
                matched_facets.append({
                    "id": fid,
                    "type": info['type'],
                    "desc": str(info)
                })
        
        print(f"  Matched Facets: {len(matched_facets)}")
        for mf in matched_facets:
            print(f"    - Facet {mf['id']} ({mf['type']})")

        results.append({
            "probe": probe['name'],
            "vanishing": vanishing,
            "diverging": diverging,
            "matched_facets": [m['id'] for m in matched_facets]
        })

    # Save
    with open(os.path.join("RESULTS", "limit_to_facet_matches_n6.json"), 'w') as f:
        json.dump(results, f, indent=2)
        
if __name__ == "__main__":
    run_probes()

