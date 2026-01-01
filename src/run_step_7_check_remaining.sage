# run_step_7_check_remaining.sage
import sys, os, time
from sage.all import *

load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_search.sage")
load("src/dcp_intersect.sage")

def main():
    print("=" * 70)
    print("STEP 7: CHECK REMAINING BOUNDARIES")
    print("=" * 70)

    # Load candidate
    input_file = "dcp_cache/boundary_12_S3xS3_s45.sobj"
    if not os.path.exists(input_file):
        log("Input file not found.")
        return
        
    Vinv_cur = load(input_file)
    log(f"Loaded candidate dim={len(Vinv_cur)}")
    
    # Load OS3
    data = load(cache_path("OS3_data.sobj"))
    C = data['C']
    triples = data['triples']
    m = len(C)
    
    # Load Flats
    Gm, full_mask = get_connected_flats(None)
    Gset = set(int(x) for x in Gm)
    incompatible = get_incompatible_pairs(None, Gm, Gset)
    triple_sign = globals().get('triple_sign') # Ensure available

    # Iterate all 2-particle channels
    # We already did (1,2) and (4,5)?
    # Let's just do ALL of them to be sure.
    
    all_2pt = []
    for ij in combinations(range(1, 7), 2):
        all_2pt.append(tuple(sorted(ij)))
        
    log(f"Checking {len(all_2pt)} 2-particle boundaries...")
    
    for S in all_2pt:
        Sidx = C.index(S)
        # S is 2-particle, Left=S, Right=Comp
        Left = set(S)
        Right = set(range(1, 7)) - Left
        
        Smask = PY1 << Sidx
        must_idx = [i for i, fm in enumerate(Gm) if int(fm) == Smask]
        Lonly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
        Ronly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)
        
        boundary_key = f"2pt_{S[0]}_{S[1]}"
        sols, best_size = run_search(Gm, Gset, full_mask, incompatible, must_idx, cache_key=boundary_key)
        
        if not sols:
            log(f"Boundary {S}: No charts found.")
            continue
            
        sol, sz = sols[0]
        res = scan_chart_exact_smallD(sol, Vinv_cur, triples, m, Lonly, Ronly, Smask, triple_sign)
        
        status = res.get('status')
        null_dim = res.get('null_dim')
        bad = res.get('bad')
        
        log(f"Boundary {S}: {status} dim={null_dim} bad={bad}")
        
        if status != 'HIT' or null_dim == 0:
            log("CANDIDATE KILLED!")
            return

    log("Candidate survived ALL 2-particle boundaries!")

if __name__ == "__main__":
    main()









