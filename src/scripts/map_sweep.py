import numpy as np
import sys
import os
import json
import itertools

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from posgeom.saddle_pushforward import compute_pushforward_saddle
from posgeom.forest_polytope import get_forest_exponents
from chy_oracle.kinematics_samples import sample_spinors_from_twistor
from posgeom.physics_map import eval_edge_vars_from_spinors
from chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian

def run_map_sweep_n6():
    print("Running Map Sweep for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Setup Polytope Data
    print("Building Forest Exponents...")
    exponents_list, edge_order = get_forest_exponents(n, roots)
    exponents = np.array(exponents_list)
    coeffs = np.ones(len(exponents))
    
    # Intrinsic Projection
    print("Projecting to Intrinsic Space...")
    v0 = exponents[0]
    diffs = exponents[1:] - v0
    
    # SVD for basis
    U, S, Vt = np.linalg.svd(diffs.T)
    rank = np.sum(S > 1e-10)
    print(f"Intrinsic Dimension: {rank}")
    
    basis = U[:, :rank]
    proj_exponents = np.dot(exponents - v0, basis)
    
    # 2. Physics Data Loop
    trials = 10
    results = []
    
    for t in range(trials):
        try:
            # Kinematics
            lambdas, tildes = sample_spinors_from_twistor(n=n, seed=200+t)
            x_ref = [1, 0]
            y_ref = [0, 1]
            
            # Edge variables z_ij
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x_ref, y_ref)
        except ValueError as e:
            print(f"Trial {t}: Kinematics error {e}")
            continue
        
        # We need to construct the "z" vector corresponding to edge_order
        # edge_order contains (u, v) tuples
        z_vec = []
        for u, v in edge_order:
            key = f"z_{u}_{v}"
            # z_map keys might be tuples or strings?
            # physics_map.py returns dict with keys (i,j) and 'z_i_j'
            if key in z_map:
                z_vec.append(float(np.abs(z_map[key]))) # Use Magnitude? Or Complex?
                # Saddle point usually for positive geometry -> real positive parameters
                # We'll take absolute value for this "positive" test.
                # If phases matter, this is a separate issue (complex saddle).
            else:
                 # Try tuple
                 if (u,v) in z_map:
                      z_vec.append(float(np.abs(z_map[(u,v)])))
                 else:
                      print(f"Missing {key}")
                      z_vec.append(1.0)
        
        z_vec = np.array(z_vec)
        
        # Candidate A: Moment Map of z
        # X_A = sum z^v * v / sum z^v
        # log_z = log(z_vec)
        # Weights ~ z^v = exp(v . log_z)
        
        # Compute weights manually
        log_z = np.log(z_vec + 1e-16)
        log_monomials = np.dot(exponents, log_z)
        shift = np.max(log_monomials)
        w = np.exp(log_monomials - shift)
        w_sum = np.sum(w)
        w_norm = w / w_sum
        
        # Target Point in Intrinsic Space
        X_target_A = np.dot(w_norm, proj_exponents)
        
        # Candidate B: Legendre Dual (Approx)
        # Assume t_vec ~ log(z_vec)
        # If z_e = e^{t_e}, then X = MomentMap(t).
        # This is effectively Candidate A.
        # "Legendre dual covector" implies we use log z as the dual coordinate directly.
        # So X_target is the moment map image of log z.
        # Which is exactly X_target_A.
        # Wait, the prompt distinguishes them?
        # "Candidate A: X = MomentMap(z)" -> This is what I computed.
        # "Candidate B: ... t is a plausible W-coordinate... X = grad Psi(t)" -> This is also MomentMap(t).
        # If t = log z, they are the same.
        # Maybe Candidate B implies optimizing t?
        # Let's stick to X_target_A for now.
        
        # Compute Pushforward Omega(X_target)
        try:
            omega_val = compute_pushforward_saddle(X_target_A, coeffs, proj_exponents)
        except Exception as e:
            print(f"Trial {t}: Saddle fail {e}")
            omega_val = np.nan
            
        # Get Gravity Amplitude
        # We need magnitude
        M_amp, _ = reconstruct_mhv_from_laplacian(lambdas, tildes, x_ref, y_ref, roots=roots)
        M_mag = float(np.abs(M_amp)) if M_amp is not None else np.nan
        
        # Forest Polynomial Value F(z)
        F_val = np.sum(w) * np.exp(shift) # sum z^v
        
        # Ratio
        ratio_omega = M_mag / omega_val if omega_val > 1e-12 else 0
        ratio_poly = M_mag / F_val if F_val > 1e-12 else 0
        
        results.append({
            "trial": t,
            "M_mag": M_mag,
            "Omega_val": omega_val,
            "F_val": F_val,
            "Ratio_Omega": ratio_omega,
            "Ratio_Poly": ratio_poly
        })
        
        print(f"Trial {t}: Ratio(Amp/Omega)={ratio_omega:.4e}, Ratio(Amp/Poly)={ratio_poly:.4e}")
        
    # Analyze
    ratios_omega = [r['Ratio_Omega'] for r in results if not np.isnan(r['Ratio_Omega'])]
    cv_omega = np.std(ratios_omega) / np.mean(ratios_omega) if ratios_omega else np.nan
    
    print(f"\nResults Summary:")
    print(f"CV of Amp/Omega: {cv_omega:.4f}")
    
    # Save
    with open('phase_u/map_sweep_results_n6.json', 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    run_map_sweep_n6()

