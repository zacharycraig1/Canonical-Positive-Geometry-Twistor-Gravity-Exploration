# run_step_4_oracle_check.sage
import sys, os, time
from sage.all import *

# Load modules
load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_invariants.sage")
load("src/dcp_search.sage")
load("src/dcp_intersect.sage")
load("src/dcp_eval.sage")  # For eval_form_on_spinors
load("src/oracle_gravity_mhv6.sage") # For oracle_M6
load("src/spinor_sampling.sage") # For sampling

def dcp_log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")

def main():
    print("=" * 70)
    print("STEP 4: ORACLE VERIFICATION")
    print("=" * 70)
    
    candidate_file = "dcp_cache/boundary_12_S3xS3_s45.sobj"
    
    if not os.path.exists(candidate_file):
        dcp_log(f"Candidate file not found: {candidate_file}")
        return

    Vinv_final = load(candidate_file)
    if len(Vinv_final) != 1:
        dcp_log(f"Candidate dim={len(Vinv_final)}, expected 1.")
        return
        
    cand_vec = Vinv_final[0]
    dcp_log(f"Loaded candidate vector (len={len(cand_vec)}).")

    # Load OS3 data to ensure basis consistency
    try:
        data = load(cache_path("OS3_data.sobj"))
        C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']
        dcp_log("Loaded OS3 data.")
    except Exception as e:
        dcp_log(f"Error loading OS3 data: {e}")
        return

    # Check sparsity again
    nz = [(i, c) for i, c in cand_vec.dict().items()]
    dcp_log(f"Sparsity: {len(nz)} non-zero entries.")

    # -------------------------------------------------------------------------
    # Verification Loop
    # -------------------------------------------------------------------------
    N_SAMPLES = 5
    ratios = []
    
    dcp_log(f"\nChecking SYMMETRIZED candidate against Oracle on {N_SAMPLES} points...")
    
    # Precompute coset representatives for S6 / (S3 x S3)
    # Maps {1,2,3} -> S, {4,5,6} -> S^c
    coset_reps = []
    for S in combinations(range(1, 7), 3):
        S = sorted(list(S))
        comp = sorted(list(set(range(1, 7)) - set(S)))
        # Map 1,2,3 -> S[0],S[1],S[2] and 4,5,6 -> comp[0],comp[1],comp[2]
        # Permutation array p such that p[i] = image of i+1
        # But Sage permutation constructor takes list of images.
        # p = [S[0], S[1], S[2], comp[0], comp[1], comp[2]]
        p = S + comp
        coset_reps.append(Permutation(p))
    dcp_log(f"Generated {len(coset_reps)} coset representatives.")

    def permute_spinors(l, lt, p):
        # p(i) is image of i (1-based).
        # We want the particle at index i to be the one that was at index p(i+1)-1?
        # No, "evaluating on permuted kinematics".
        # f(sigma * lambda).
        # Usually means particle 1 is now what particle sigma(1) was.
        # So new_l[0] = l[p(1)-1].
        return ([l[p(i+1)-1] for i in range(6)], 
                [lt[p(i+1)-1] for i in range(6)])

    consistent_count = 0
    
    for i in range(N_SAMPLES):
        seed = 1000 + i
        sample = sample_spinor_helicity_conserving(n=6, seed=seed)
        if sample is None: continue
        lambdas, tilde_lambdas = sample
        
        # Symmetrized Evaluation
        val_sym = QQ(0)
        for p in coset_reps:
            l_p, lt_p = permute_spinors(lambdas, tilde_lambdas, p)
            val = eval_form_on_spinors(cand_vec, l_p, lt_p, C, triples)
            val_sym += val
            
        # Oracle
        res_oracle = oracle_M6(lambdas, tilde_lambdas, convention="hodges_reduced")
        val_oracle = res_oracle.get("value")
        
        if val_oracle is None or val_oracle == 0:
            dcp_log(f"  Sample {i}: Oracle invalid")
            continue
            
        ratio = val_sym / val_oracle
        ratios.append(ratio)
        
        dcp_log(f"  Sample {i}: sym_cand={val_sym}, oracle={val_oracle}")
        dcp_log(f"    Ratio={ratio}")

    # Summary
    print("-" * 30)
    if not ratios: return

    first = ratios[0]
    if all(r == first for r in ratios):
        dcp_log(f"SUCCESS! Symmetrized candidate matches Oracle! Ratio = {first}")
    else:
        dcp_log("FAILURE. Ratios vary.")
        dcp_log(f"Ratios: {ratios}")

if __name__ == "__main__":
    main()

