
# src/run/run_true_pullback_scan_S123.sage
import sys
import os
sys.path.append(os.getcwd())

# Load modules
load('src/dcp_common.sage')
load('src/dcp/connected_flats.sage')
load('src/os/s6_invariants.sage')
load('src/dcp/nested_set_sampler_seeded.sage')
load('src/dcp/chart_pullback_true.sage')

def log(msg):
    print(f"[TRUE_DCP] {msg}")

def run_scan():
    log("=======================================================")
    log("TRUE DCP PULLBACK SCAN")
    log("=======================================================")

    # 1. Setup Environment
    C, M6 = build_M6_matroid()
    Gm, full_mask = get_connected_flats(M6)
    incompatible = get_incompatible_pairs(M6, Gm, set(Gm))
    triples, Vbasis = build_OS3_data(C, M6)
    
    # 2. Load Invariants
    Vinv_vectors = compute_S6_invariants(C, triples, Vbasis)
    
    # 3. Define Seed and Boundary
    chan_to_idx = {ch: i for i, ch in enumerate(C)}
    
    def get_flat_mask_for_particles(particles):
        p_list = sorted(list(particles))
        if len(p_list) == 2: ch = tuple(p_list)
        elif len(p_list) == 3: ch = canonical_three_subset(set(p_list), 6)
        else: return 0
        if ch not in chan_to_idx: return 0
        idx = chan_to_idx[ch]
        return set_to_mask(M6.closure([idx]))

    seed_tuples = [(1,2), (1,2,3), (5,6)]
    seed_masks = []
    for tup in seed_tuples:
        m = get_flat_mask_for_particles(tup)
        seed_masks.append(m)
            
    S_tuple = (1,2,3)
    S_mask = get_flat_mask_for_particles(S_tuple)
    
    L_set = {1,2,3}
    R_set = {4,5,6}
    Lonly, Ronly = build_LR_masks(C, L_set, R_set)
    
    # 4. Test on Seed First
    log("Testing SEED ONLY (Partial Chart)...")
    # Using true pullback on the seed
    # Note: rank of seed is 3. Adapted basis will complete it to 10.
    
    try:
        result_seed = compute_dcp_residue(seed_masks, Vinv_vectors, triples, M6, C, S_mask)
        log(f"Seed Result: {result_seed.get('status')} Bad={result_seed.get('bad_count')}")
    except Exception as e:
        log(f"Seed Test Failed with error: {e}")
        import traceback
        traceback.print_exc()
        result_seed = None

    # 5. Test Completions
    log("Generating completions...")
    samples = sample_nested_sets_seeded(Gm, incompatible, seed_masks, n_samples=10, rng_seed=42)
    log(f"Generated {len(samples)} samples.")
    
    hits = []
    
    for i, N in enumerate(samples):
        log(f"Sample {i}: size={len(N)}")
        try:
            result = compute_dcp_residue(N, Vinv_vectors, triples, M6, C, S_mask)
            status = result.get('status')
            bad = result.get('bad_count')
            log(f"  Status: {status}, Bad: {bad}")
            
            if status == 'HIT':
                hits.append((i, result))
                log("  HIT FOUND!")
        except Exception as e:
            log(f"  Error on sample {i}: {e}")
            
    # 6. Report
    log("\nSUMMARY:")
    if result_seed:
        log(f"Seed: {result_seed.get('status')}")
    log(f"Hits in completions: {len(hits)}")
    
    # Save a hit if found
    if hits:
        hit = hits[0][1]
        save(hit, os.path.join(CACHE_DIR, "TRUE_PULLBACK_HIT.sobj"))
    elif result_seed and result_seed.get('status') == 'HIT':
        save(result_seed, os.path.join(CACHE_DIR, "TRUE_PULLBACK_HIT_SEED.sobj"))

if __name__ == '__main__':
    run_scan()






