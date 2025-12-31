#!/usr/bin/env sage
from sage.all import *
import os
import time

# Load libraries
load('src/hodges.sage')
load('src/klt.sage')

def generate_data_batch(batch_id, batch_size=100):
    print(f"Generating Data Batch {batch_id} (Size: {batch_size})")
    
    # Load basis metadata
    if not os.path.exists("klt_search/basis_metadata.sobj"):
        print("Error: Metadata not found. Run 01_setup_basis.sage first.")
        return
        
    basis_data = load("klt_search/basis_metadata.sobj")
    perms = basis_data['perms']
    monomials = basis_data['monomials']
    # basis_names = basis_data['basis_names'] # Not strictly needed if we assume index order
    
    # Pre-compute fixed legs for Parke-Taylor
    # Fixed: 1, 5, 6 (indices 0, 4, 5)
    # Permuted: 2, 3, 4 (indices 1, 2, 3)
    
    data_points = []
    
    count = 0
    attempts = 0
    start_time = time.time()
    
    while count < batch_size:
        attempts += 1
        
        # 1. Generate Random Kinematics
        # using check_domain=False for speed, but need to be careful
        twistor = MomentumTwistor(n=6, check_domain=False)
        
        try:
            # 2. Compute Oracle Amplitude (Hodges Reduced)
            # We must use the reduced form and normalize it to match KLT scaling.
            # KLT scales as t^-2. Hodges Reduced scales as t^-18.
            # Normalization factor: - <0 1>^8 (scales as t^16).
            
            # Use Adapter for reduced function
            load('src/kinematics_map.sage')
            load('src/spinor_helicity.sage')
            
            lambdas, tilde_lambdas, _ = extract_spinors_from_twistor(twistor)
            if lambdas is None or tilde_lambdas is None: continue
            
            adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
            
            # hodges_6pt_mhv_reduced is loaded by src/hodges.sage
            
            amp_hodges_red, reason = hodges_6pt_mhv_reduced(adapter)
            if amp_hodges_red is None: continue
            
            ang01 = twistor.get_angle(0, 1)
            if ang01 == 0: continue
            
            norm_factor = - (ang01 ** 8)
            target_val = amp_hodges_red * norm_factor
            
            # 3. Compute Basis Variables (Mandelstams)
            # Order: s12, s23, s34, s45, s56, s61, s123, s234, s345
            s_vars = []
            
            # 2-pt
            pairs = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
            for i, j in pairs:
                val = mandelstam_invariant(twistor, i, j)
                if val is None: raise ValueError("Bad s_ij")
                s_vars.append(val)
                
            # 3-pt
            triples = [(0,1,2), (1,2,3), (2,3,4)]
            for i, j, k in triples:
                s_ij = mandelstam_invariant(twistor, i, j)
                s_jk = mandelstam_invariant(twistor, j, k)
                s_ki = mandelstam_invariant(twistor, k, i)
                if s_ij is None: raise ValueError("Bad s_ijk")
                s_vars.append(s_ij + s_jk + s_ki)
                
            # 4. Compute Monomial Values
            # Monomials are tuples of exponents corresponding to s_vars indices
            mon_vals = []
            for exps in monomials:
                val = 1
                for idx, e in enumerate(exps):
                    if e > 0:
                        val *= (s_vars[idx] ** e)
                mon_vals.append(val)
            
            # 5. Compute Parke-Taylor Basis Values
            pt_vals = []
            for p in perms:
                # Ordering: 1, p..., 5, 6 -> [0] + p + [4, 5]
                order = [0] + list(p) + [4, 5]
                val = parke_taylor_6pt_mhv(twistor, order)
                if val is None: raise ValueError("Bad PT")
                pt_vals.append(val)
                
            # Store the point data
            # We don't build the full row here to save space/time, just store the components.
            # Row construction happens in solver.
            # Structure: { 'target': val, 'mon_vals': [...], 'pt_vals': [...] }
            
            point_data = {
                'target': target_val,
                'mon_vals': mon_vals,
                'pt_vals': pt_vals
            }
            data_points.append(point_data)
            count += 1
            
        except (ValueError, ArithmeticError):
            continue
            
    end_time = time.time()
    print(f"Generated {count} points in {end_time - start_time:.2f}s")
    
    # Save batch
    save(data_points, f"klt_search/data_batch_{batch_id}.sobj")
    
if __name__ == "__main__":
    import sys
    # Usage: sage 02_generate_data.sage [batch_id] [batch_size]
    args = sys.argv[1:]
    batch_id = int(args[0]) if len(args) > 0 else 0
    batch_size = int(args[1]) if len(args) > 1 else 100
    
    generate_data_batch(batch_id, batch_size)

