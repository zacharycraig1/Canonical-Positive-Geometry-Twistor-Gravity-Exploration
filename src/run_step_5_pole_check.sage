# run_step_5_pole_check.sage
import sys, os
from sage.all import *

load("src/dcp_common.sage")
load("src/dcp_eval.sage")
load("src/oracle_gravity_mhv6.sage")
load("src/boundary_sampler.sage")

def main():
    print("=" * 70)
    print("STEP 5: POLE BEHAVIOR CHECK")
    print("=" * 70)

    candidate_file = "dcp_cache/boundary_12_S3xS3_s45.sobj"
    if not os.path.exists(candidate_file):
        print(f"Candidate file not found: {candidate_file}")
        return
    Vinv_final = load(candidate_file)
    cand_vec = Vinv_final[0]
    
    data = load(cache_path("OS3_data.sobj"))
    C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']

    # Define channel to probe
    # s_12 (indices 0 and 1) -> channel (1, 2)
    # In C, (1,2) is typically early.
    target_channel = (1, 2)
    
    print(f"\nProbing pole at s_{target_channel} -> 0")
    
    # Generate a point near the boundary
    # sample_near_boundary(n, s_channel_indices, epsilon)
    # But boundary_sampler usually works with indices 0..n-1
    # Channel (1,2) corresponds to indices {0, 1}
    
    # We will manually shift using boundary_sampler functions
    # 1. Start with random point
    res = sample_spinor_helicity_conserving(n=6, seed=42)
    if not res: return
    lambdas, tildes = res
    
    # 2. Shift to make s_12 small
    # Shift params: probe channel S={0,1}, shift indices (i,j)=(0,3) to match test
    S_indices = [0, 1]
    shift_indices = (0, 3) 
    
    try:
        z_pole = solve_boundary_shift(lambdas, tildes, S_indices, shift_indices)
    except Exception as e:
        print(f"Shift failed: {e}")
        return
        
    print(f"Found pole at z = {z_pole}")
    
    epsilons = [1/10, 1/100, 1/1000, 1/10000]
    
    for eps in epsilons:
        z = z_pole + eps
        l_new, lt_new = apply_shift(lambdas, tildes, shift_indices[0], shift_indices[1], z)
        
        # Check s_12
        s = spinors_to_sij(l_new, lt_new)
        s12 = s[(0,1)]
        # print(f"  eps={float(eps):.1e}, s_12={complex(s12):.1e}")
        
        # Eval
        val_cand = eval_form_on_spinors(cand_vec, l_new, lt_new, C, triples)
        
        res_oracle = oracle_M6(l_new, lt_new, convention="hodges_reduced")
        val_oracle = res_oracle.get("value")
        
        if val_oracle is None:
            print("  Oracle failed")
            continue
            
        ratio = val_cand / val_oracle
        
        # Check scaling
        # val * s12 should be constant near pole
        residue_cand = val_cand * s12
        residue_oracle = val_oracle * s12
        
        print(f"  eps={float(eps):.5f} | s12={float(s12):.5f}")
        print(f"    Cand Val: {float(abs(val_cand)):.2e}  (Res: {float(abs(residue_cand)):.4f})")
        print(f"    Orac Val: {float(abs(val_oracle)):.2e}  (Res: {float(abs(residue_oracle)):.4f})")
        print(f"    Ratio: {float(abs(ratio)):.4f}")

if __name__ == "__main__":
    main()

