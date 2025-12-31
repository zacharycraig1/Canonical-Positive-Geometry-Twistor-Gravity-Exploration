# proof_residue_audit.sage
import sys, os, time, json
from sage.all import *

load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_search.sage")
load("src/dcp_intersect.sage")

def json_converter(obj):
    if hasattr(obj, 'is_integer') and obj.is_integer(): return int(obj)
    if hasattr(obj, 'is_real') and obj.is_real(): return float(obj)
    if isinstance(obj, (int, Integer)): return int(obj)
    if isinstance(obj, (float)): return float(obj)
    if isinstance(obj, list): return [json_converter(x) for x in obj]
    if isinstance(obj, dict): return {str(k): json_converter(v) for k,v in obj.items()}
    return str(obj)

def main():
    print("=" * 70)
    print("PROOF A1: RESIDUE OPERATOR AUDIT")
    print("=" * 70)
    
    # 1. Load OS3 and Flats
    try:
        data = load(cache_path("OS3_data.sobj"))
        C = data['C']
        m = len(C)
        triples = data['triples']
        Gm, full_mask = get_connected_flats(None)
        Gset = set(int(x) for x in Gm)
        incompatible = get_incompatible_pairs(None, Gm, Gset)
        sign_table = get_sign_table()
        triple_sign = globals().get('triple_sign')
    except Exception as e:
        log(f"Error loading prerequisites: {e}")
        return

    # 2. Pick divisor s_123
    S = (1, 2, 3)
    S_can = canonical_three_subset(S, 6)
    S_ch = channel_of_three_subset(S_can, 6)
    Sidx = C.index(S_ch)
    Smask = PY1 << Sidx
    
    must_idx = [i for i, fm in enumerate(Gm) if int(fm) == Smask]
    
    # Left/Right for factorization (needed for scan logic)
    Left = set(S_can)
    Right = set(range(1, 7)) - Left
    Lonly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
    Ronly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)

    # 3. Find 2 different charts containing S
    # We run search with random seed until we get enough unique solutions
    # that differ significantly.
    
    # Just run search once, it usually returns 100+ charts.
    log(f"Searching for charts containing {S}...")
    
    solutions, best_size = run_search(Gm, Gset, full_mask, incompatible, must_idx)
    
    if len(solutions) < 2:
        log("Not enough charts found to compare!")
        return
        
    chart1, sz1 = solutions[0]
    chart2, sz2 = solutions[-1] # Pick last to be different
    
    if chart1 == chart2:
        log("Charts identical?")
        if len(solutions) > 1:
            chart2, sz2 = solutions[1]
    
    log(f"Chart 1 (size {sz1})")
    log(f"Chart 2 (size {sz2})")
    
    # 4. Generate constraint matrices for both
    # We use a dummy Vinv (identity) to get the full constraint matrix on coefficients
    # d = global dimension (2008)
    
    # This might be too large for full constraint matrix (2008 cols).
    # Instead, let's verify on a random subspace or just check dim(Null).
    
    # Let's take Vinv = Identity(2008) is too big.
    # Let's take Vinv = S3xS3 basis (58 dims) or random 50 vectors?
    # Better: S3xS3 basis, since that's where we work.
    
    try:
        V_basis = load("dcp_cache/basis_S3xS3.sobj")
        log(f"Using S3xS3 basis (dim {len(V_basis)}) for audit.")
    except:
        log("S3xS3 basis not found, creating random...")
        V_basis = [vector(QQ, [1 if i==j else 0 for i in range(len(triples))]) for j in range(10)]
        
    # Get constraints for Chart 1
    res1 = scan_chart_exact_smallD(chart1, V_basis, triples, m, Lonly, Ronly, Smask, triple_sign, return_basis=True)
    null1 = res1.get('null_dim')
    bad1 = res1.get('bad')
    log(f"Chart 1: dim(Null) = {null1}, bad pairs = {bad1}")
    
    # Get constraints for Chart 2
    res2 = scan_chart_exact_smallD(chart2, V_basis, triples, m, Lonly, Ronly, Smask, triple_sign, return_basis=True)
    null2 = res2.get('null_dim')
    bad2 = res2.get('bad')
    log(f"Chart 2: dim(Null) = {null2}, bad pairs = {bad2}")
    
    # 5. Compare Nullspaces
    # The nullspaces should be identical (subspace of V_basis).
    # Nullspace is represented by basis V_new = V_basis * Cmat.
    # Cmat1 = from res1['_B']
    # Cmat2 = from res2['_B']
    
    B1 = res1.get('_B')
    sc1 = res1.get('_scales')
    C1 = _coeff_matrix_from_scan_B(B1, sc1) # d x k1
    
    B2 = res2.get('_B')
    sc2 = res2.get('_scales')
    C2 = _coeff_matrix_from_scan_B(B2, sc2) # d x k2
    
    # Check if col_space(C1) == col_space(C2)
    # They are subspaces of Q^d (d=58).
    
    if null1 != null2:
        log("FAIL: Dimensions differ!")
    else:
        # Check equality
        # Rank of [C1 | C2] should be null1
        M_aug = block_matrix([[C1, C2]])
        rank_aug = M_aug.rank()
        if rank_aug == null1:
            log("SUCCESS: Nullspaces are identical.")
            status = "PASS"
        else:
            log(f"FAIL: Nullspaces differ. Joint rank = {rank_aug} (expected {null1})")
            status = "FAIL"

    # Export
    out = {
        "divisor": str(S),
        "chart1_size": int(sz1),
        "chart2_size": int(sz2),
        "dim_subspace": len(V_basis),
        "null_dim_chart1": int(null1),
        "null_dim_chart2": int(null2),
        "status": status
    }
    
    with open("proof/residue_audit_s123.json", "w") as f:
        json.dump(json_converter(out), f, indent=2)
        
    with open("proof/residue_audit_s123.txt", "w") as f:
        f.write(f"Residue Audit for {S}\n")
        f.write(f"Chart 1 Size: {sz1}\n")
        f.write(f"Chart 2 Size: {sz2}\n")
        f.write(f"Subspace Dim: {len(V_basis)}\n")
        f.write(f"Chart 1 Null Dim: {null1}\n")
        f.write(f"Chart 2 Null Dim: {null2}\n")
        f.write(f"Status: {status}\n")

if __name__ == "__main__":
    main()

