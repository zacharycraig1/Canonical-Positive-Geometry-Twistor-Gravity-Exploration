# run_step_3_verify.sage
import sys, os, time
from sage.all import *

# Load modules
load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_invariants.sage")
load("src/dcp_search.sage") # for build_candidate_vec
load("src/dcp_intersect.sage") # for helpers

def main():
    args = sys.argv[1:]
    if len(args) < 1:
        # Default to the one we just found
        candidate_file = "dcp_cache/boundary_12_S3xS3_s45.sobj"
    else:
        candidate_file = args[0]

    print("=" * 70)
    print(f"STEP 3: VERIFY CANDIDATE")
    print(f"Input: {candidate_file}")
    print("=" * 70)

    if not os.path.exists(candidate_file):
        log(f"Candidate file not found: {candidate_file}")
        return

    Vinv_final = load(candidate_file)
    if len(Vinv_final) != 1:
        log(f"Candidate has dimension {len(Vinv_final)}, expected 1.")
        return

    # Load OS3 data
    try:
        data = load(cache_path("OS3_data.sobj"))
        C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']
    except Exception as e:
        log(f"Error loading OS3 data: {e}")
        return

    # For verification we need Tmat, but we ran this step-by-step so Tmat is lost?
    # Actually verify_dim1_candidate takes alpha0, Vinv_final, Tmat...
    # But Vinv_final[0] IS the candidate vector in OS3 space (if we treated it as such).
    # Wait, Vinv_final contains vectors in OS3 space (combined from previous steps).
    # dcp_intersect.sage: Vinv_cur = _combine_invariant_basis(Vinv_cur, Cmat) returns list of vectors in OS3 basis.
    # So Vinv_final[0] is the vector.
    
    cand_vec = Vinv_final[0]
    
    # Check sparsity (number of triples)
    nz = [(i, c) for i, c in cand_vec.dict().items()]
    log(f"Candidate sparsity: {len(nz)} non-zero triples")
    
    # Check S6 invariance
    log("Checking S6 invariance...")
    # Load S6 invariants to project against? Or just generate them.
    # We can use check_s6_invariance helper if we wrote one.
    # Or just project onto S6 space.
    
    # Compute S6 invariants again (or load)
    try:
        V_s6 = load(cache_path("basis_S6.sobj"))
    except:
        V_s6 = compute_S6_invariants(C, triples, Vbasis)
        
    # Check if cand_vec is in span(V_s6)
    # M = matrix([cand_vec])
    # B = matrix(V_s6)
    # Rank check
    M_s6 = matrix(QQ, V_s6)
    M_aug = M_s6.stack(vector(QQ, cand_vec))
    if M_aug.rank() == M_s6.rank():
        log("SUCCESS: Candidate IS S6 invariant!")
    else:
        log("FAILURE: Candidate is NOT S6 invariant.")
        
    # Check against Oracle (Optional/Advanced)
    # We need to import oracle and evaluate.
    # Not strictly required for this step, but good confirmation.
    
    log("Verification complete.")

if __name__ == "__main__":
    main()










