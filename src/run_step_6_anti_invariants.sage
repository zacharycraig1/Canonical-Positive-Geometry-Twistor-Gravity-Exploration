# run_step_6_anti_invariants.sage
import sys, os, time
from sage.all import *

load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_search.sage")
load("src/dcp_intersect.sage")

def compute_anti_invariants_from_generators(C, triples, Vbasis, generators):
    """Compute ANTI-invariants (sign rep) under a subgroup."""
    dimV, Wdim = Vbasis.nrows(), Vbasis.ncols()

    def sign_sort(triple):
        arr = list(triple)
        sign = 1
        for i in range(2):
            for j in range(2-i):
                if arr[j] > arr[j+1]:
                    arr[j], arr[j+1] = arr[j+1], arr[j]
                    sign = -sign
        return sign, tuple(arr)

    B = matrix(QQ, Wdim, dimV, sparse=True)
    for j in range(dimV):
        for i, val in vector(QQ, Vbasis.row(j)).dict().items():
            B[i, j] = val

    chan_index = {C[i]: i for i in range(len(C))}
    triple_to_row = {triples[i]: i for i in range(len(triples))}

    for t, perm in enumerate(generators):
        img_row, img_sgn = [0]*len(triples), [0]*len(triples)
        for i, tri in enumerate(triples):
            imgs = [chan_index[_perm_ch(C[idx], perm)] for idx in tri]
            if len(set(imgs)) < 3:
                img_row[i], img_sgn[i] = i, 0
            else:
                sgn, tri2 = sign_sort(tuple(imgs))
                img_row[i], img_sgn[i] = triple_to_row[tri2], sgn

        PB = matrix(QQ, Wdim, B.ncols(), sparse=True)
        for (r, c), val in B.dict().items():
            if img_sgn[r]:
                PB[img_row[r], c] += img_sgn[r] * val

        # For anti-invariants: PB * v = sgn(perm) * v
        # Generators are transpositions (odd), so sgn(perm) = -1
        # We want PB * v = -v => (PB + B) * v = 0
        
        perm_sign = perm.sign()
        if perm_sign == -1:
            diff = PB + B
        else:
            diff = PB - B # Should be invariant if even permutation? 
            # Usually we use generators that are transpositions.
            
        kernel_basis = fast_right_kernel(diff)
        B = B * kernel_basis.transpose()
        log(f"  anti-inv gen {t+1}/{len(generators)}: dim = {B.ncols()}")

    return [B.column(i) for i in range(B.ncols())]

def main():
    print("=" * 70)
    print("STEP 6: ANTI-INVARIANT SEARCH")
    print("=" * 70)

    try:
        data = load(cache_path("OS3_data.sobj"))
        C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']
        m = len(C)
        triples_array = np.array(triples, dtype=np.int32)
    except Exception as e:
        log(f"Error loading OS3 data: {e}")
        return

    # S6 Anti-Invariants
    log("Computing S6 ANTI-invariants...")
    gens = [SymmetricGroup(6)((i, i+1)) for i in range(1, 6)]
    V_s6_anti = compute_anti_invariants_from_generators(C, triples, Vbasis, gens)
    log(f"  S6 Anti-Invariants: {len(V_s6_anti)} vectors")
    
    if len(V_s6_anti) == 0:
        log("No S6 anti-invariants found. Stopping.")
        return

    # If found, try intersection
    Gm, full_mask = get_connected_flats(None)
    Gset = set(int(x) for x in Gm)
    incompatible = get_incompatible_pairs(None, Gm, Gset)
    sign_table = get_sign_table() # Need to init
    triple_sign = globals()['triple_sign']

    # Try boundary (1,2,3)
    log("\nIntersecting S6-Anti with boundary (1,2,3)...")
    # ... logic from run_step_2 ...
    S = (1, 2, 3)
    # ... setup ...
    S_can = canonical_three_subset(S, 6)
    S_ch = channel_of_three_subset(S_can, 6)
    Sidx = C.index(S_ch)
    Left = set(S_can); Right = set(range(1,7)) - Left
    Smask = PY1 << Sidx
    must_idx = [i for i, fm in enumerate(Gm) if int(fm) == Smask]
    Lonly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
    Ronly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)
    
    boundary_key = "_".join(map(str, sorted(list(S))))
    sols, best_size = run_search(Gm, Gset, full_mask, incompatible, must_idx, cache_key=boundary_key)
    
    if sols:
        sol, sz = sols[0]
        res = scan_chart_exact_smallD(sol, V_s6_anti, triples, m, Lonly, Ronly, Smask, triple_sign, return_basis=True)
        log(f"Result: {res.get('status')} dim={res.get('null_dim')}")
        
        if res.get('null_dim') > 0:
            log("FOUND SURVIVING S6 ANTI-INVARIANTS!")
            # Save candidate
            B_mat = res.get('_B')
            scales = res.get('_scales')
            Cmat = _coeff_matrix_from_scan_B(B_mat, scales)
            V_survived = _combine_invariant_basis(V_s6_anti, Cmat)
            save(V_survived, cache_path("candidate_S6_anti.sobj"))
            log("Saved to candidate_S6_anti.sobj")
    else:
        log("No charts.")

if __name__ == "__main__":
    main()






