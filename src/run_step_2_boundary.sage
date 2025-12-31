# run_step_2_boundary.sage
import sys, os, time, argparse
from sage.all import *

# Load modules
load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_search.sage")
load("src/dcp_intersect.sage")

def parse_boundary(s):
    # expect "1,2,3"
    return tuple(map(int, s.split(',')))

def main():
    # Manual argument parsing since we are running inside sage script
    # Expected usage: sage run_step_2_boundary.sage BOUNDARY INPUT_BASIS OUTPUT_BASIS
    
    args = sys.argv[1:]
    # Filter out sage specific args if any (usually not passed to script logic directly but be careful)
    
    if len(args) < 3:
        print("Usage: sage run_step_2_boundary.sage 1,2,3 input_basis.sobj output_basis.sobj")
        # Try to guess or default for testing
        return

    boundary_str = args[0]
    input_file = args[1]
    output_file = args[2]

    S = parse_boundary(boundary_str)
    
    print("=" * 70)
    print(f"STEP 2: INTERSECT BOUNDARY {S}")
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print("=" * 70)

    # Load OS3 data
    try:
        data = load(cache_path("OS3_data.sobj"))
        C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']
        m = len(C)
        triples_array = np.array(triples, dtype=np.int32)
    except Exception as e:
        log(f"Error loading OS3 data: {e}. Run step 0 first.")
        return

    # Load Flats/Incompat
    Gm, full_mask = get_connected_flats(None) # Pass None assuming cache exists
    if Gm is None:
        # If cache missing, need M6
        _, M6 = build_M6_matroid()
        Gm, full_mask = get_connected_flats(M6)
        
    Gset = set(int(x) for x in Gm)
    incompatible = get_incompatible_pairs(None, Gm, Gset)

    # Load Input Basis
    try:
        Vinv_cur = load(input_file)
        log(f"Loaded {len(Vinv_cur)} vectors from {input_file}")
    except Exception as e:
        log(f"Error loading input basis: {e}")
        return

    if len(Vinv_cur) == 0:
        log("Input basis is empty. Writing empty output.")
        save([], output_file)
        return

    # Prepare for single boundary intersection
    # Logic adapted from dcp_intersect.sage
    
    S_tuple = tuple(S)
    if len(S_tuple) == 2:
        S_ch = tuple(sorted(S_tuple))
        Sidx = C.index(S_ch)
        Left = set(S_tuple)
        Right = set(range(1,7)) - Left
    elif len(S_tuple) == 3:
        S_can = canonical_three_subset(S_tuple, 6)
        S_ch = channel_of_three_subset(S_can, 6)
        Sidx = C.index(S_ch)
        Left = set(S_can)
        Right = set(range(1,7)) - Left
    else:
        log(f"Invalid boundary size: {len(S)}")
        return

    Smask = PY1 << Sidx
    must_idx = [i for i, fm in enumerate(Gm) if int(fm) == Smask]
    Lonly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
    Ronly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)

    log(f"Boundary masks: Sidx={Sidx}, |Left|={len(Left)}, |Right|={len(Right)}")

    boundary_key = "_".join(map(str, sorted(list(S))))
    sols, best_size = run_search(Gm, Gset, full_mask, incompatible, must_idx, cache_key=boundary_key)
    
    if not sols:
        log("No charts found. Empty output.")
        save([], output_file)
        return

    sol, sz = sols[0]
    res = scan_chart_exact_smallD(sol, Vinv_cur, triples, m, Lonly, Ronly, Smask, triple_sign, return_basis=True)
    
    status = res.get('status')
    null_dim = int(res.get('null_dim', 0))
    bad = res.get('bad')
    
    log(f"Result: status={status}, null_dim={null_dim}, bad={bad}")
    
    if status == 'HIT' and null_dim > 0:
        B = res.get('_B')
        scales = res.get('_scales')
        Cmat = _coeff_matrix_from_scan_B(B, scales)
        Vinv_new = _combine_invariant_basis(Vinv_cur, Cmat)
        save(Vinv_new, output_file)
        log(f"Saved {len(Vinv_new)} vectors to {output_file}")
    else:
        log("Intersection empty.")
        save([], output_file)

if __name__ == "__main__":
    main()






