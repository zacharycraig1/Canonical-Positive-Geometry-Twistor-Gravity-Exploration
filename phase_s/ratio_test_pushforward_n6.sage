
import sys
import os
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.canonical_form import CanonicalFormEvaluator
from src.posgeom.physics_map import eval_edge_vars_from_spinors

# Load Hodges modules
try:
    load('src/hodges_sigma.sage')
    load('src/hodges.sage')
    print("Loaded src/hodges.sage")
except Exception as e:
    print(f"Error loading src/hodges.sage: {e}")
    # If load fails, we might be in python mode not sage mode? 
    # But we are running with sage.ps1
    pass

def run_ratio_test():
    print("S5: Pushforward Ratio Test (n=6)...")
    
    # Check if MomentumTwistor is available
    if 'MomentumTwistor' not in globals():
        print("MomentumTwistor class not found. Import failed.")
        return

    # 2. Setup Polytope (n=6, roots=[0,1,2])
    n = 6
    roots = [0, 1, 2]
    print(f"Building Forest Polytope (n={n}, roots={roots})...")
    exponents, edge_order = get_forest_exponents(n, roots)
    P = Polyhedron(vertices=exponents)
    print(f"Polytope dim: {P.dim()}")
    
    # 3. Loop trials
    num_trials = 5
    ratios = []
    
    for trial in range(num_trials):
        seed_val = 100 + trial
        set_random_seed(seed_val)
        
        # A. Generate Kinematics
        # We need MomentumTwistor
        # Random integers for twistors to stay in QQ
        data = [[randint(-10, 10) for _ in range(4)] for _ in range(n)]
        # Fix: Pass Z explicitly
        twistor = MomentumTwistor(n=n, Z=data)
        
        # Check non-degenerate
        if any(twistor.get_angle(i, (i+1)%n) == 0 for i in range(n)):
            continue
            
        # B. Compute Amplitude
        # Need reference spinors for Hodges
        ref_spinors = sample_reference_spinors(twistor)
        if ref_spinors[0] is None:
            print("Failed to sample ref spinors")
            continue
            
        amp, status = hodges_6pt_mhv_reduced(twistor, ref_spinors)
        if amp is None:
            print(f"Hodges failed: {status}")
            continue
            
        # C. Compute Canonical Form
        # Get lambdas/tildes
        lambdas = {i: vector(QQ, twistor.get_lambda(i)) for i in range(n)}
        tildes = {i: vector(QQ, twistor.get_tilde_lambda(i)) for i in range(n)}
        
        # Physics Map needs reference spinors x, y
        # We can use the same ones used for Hodges!
        lam_x, lam_y = ref_spinors
        
        try:
            z_vals = eval_edge_vars_from_spinors(lambdas, tildes, lam_x, lam_y)
        except ValueError as e:
            print(f"Map failed: {e}")
            continue
            
        # Construct W
        W_list = [1]
        for u, v in edge_order:
            if u > v: u, v = v, u
            val = z_vals.get((u, v))
            if val is None: val = z_vals.get(f"z_{u}_{v}")
            W_list.append(val)
        W = vector(QQ, W_list)
        
        # Eval Polytope
        try:
            omega = CanonicalFormEvaluator.eval_polytope(P, W)
        except Exception as e:
            print(f"Eval Polytope failed: {e}")
            continue
            
        # D. Ratio
        if omega == 0:
            print("Omega is 0!")
            continue
            
        ratio = amp / omega
        ratios.append(ratio)
        print(f"Trial {trial}: Amp={float(amp):.4e}, Omega={float(omega):.4e}, Ratio={float(ratio):.4e}")
        
    # Check consistency
    if len(ratios) > 1:
        first = ratios[0]
        consistent = all(abs(r - first)/abs(first) < 1e-5 for r in ratios)
        if consistent:
            print(f"\n[PASS] Ratio is constant: {first}")
        else:
            print("\n[FAIL] Ratio varies.")
            print(ratios)
    else:
        print("Not enough data.")

if __name__ == "__main__":
    run_ratio_test()

