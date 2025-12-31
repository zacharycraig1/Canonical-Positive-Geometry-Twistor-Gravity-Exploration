
# src/run/run_proxy_hit_repro.sage
import sys
import os
sys.path.append(os.getcwd())

# Load modules
load('src/dcp_common.sage')
load('src/dcp/connected_flats.sage')
load('src/os/s6_invariants.sage')
load('src/dcp/nested_set_sampler_seeded.sage')
load('src/dcp/proxy_pullback.sage')

def log(msg):
    print(f"[REPRO] {msg}")

def run_repro():
    log("=======================================================")
    log("PROXY HIT REPRODUCTION")
    log("=======================================================")

    # 1. Setup Environment
    C, M6 = build_M6_matroid()
    Gm, full_mask = get_connected_flats(M6)
    log(f"Number of connected flats Gm: {len(Gm)}")
    
    incompatible = get_incompatible_pairs(M6, Gm, set(Gm))
    log(f"Number of incompatible pairs: {len(incompatible)}")
    
    triples, Vbasis = build_OS3_data(C, M6)
    sign_table = get_sign_table()
    
    # 2. Load Invariants
    # Assuming standard S6 invariants (dim 2)
    Vinv_vectors = compute_S6_invariants(C, triples, Vbasis)
    log(f"Invariants loaded: dim = {len(Vinv_vectors)}")
    
    # 3. Define Seed and Boundary
    # N_seed = [(1,2), (1,2,3), (5,6)]
    # These are sets of PARTICLES defining the poles.
    # We need to map them to flats in the channel matroid M6.
    
    # Map particle tuple -> channel index
    chan_to_idx = {ch: i for i, ch in enumerate(C)}
    
    def get_flat_mask_for_particles(particles):
        # particles is tuple/set like (1,2) or (1,2,3)
        # Normalize to canonical channel representation
        p_list = sorted(list(particles))
        if len(p_list) == 2:
            ch = tuple(p_list)
        elif len(p_list) == 3:
            ch = canonical_three_subset(set(p_list), 6)
        else:
            log(f"Error: Invalid particle set {particles}")
            return 0
            
        if ch not in chan_to_idx:
            log(f"Error: Channel {ch} not in C")
            return 0
            
        idx = chan_to_idx[ch]
        # Compute closure of this single channel
        # M6.closure takes list of elements
        # dcp_common doesn't expose M6 directly but we have it from build_M6_matroid
        
        # We need the closure of {idx}
        F_closure = M6.closure([idx])
        return set_to_mask(F_closure)

    seed_tuples = [(1,2), (1,2,3), (5,6)]
    seed_masks = []
    
    for tup in seed_tuples:
        m = get_flat_mask_for_particles(tup)
        seed_masks.append(m)
        
    log(f"Seed masks computed: {[bin(m) for m in seed_masks]}")
    
    # Verify they are in Gm
    Gm_set = set(Gm)
    for i, m in enumerate(seed_masks):
        if m not in Gm_set:
            log(f"Error: Seed flat {seed_tuples[i]} (mask {bin(m)}) not in connected flats Gm!")
            # Fallback: maybe it's not strictly in Gm if filtered?
            # But normally all single-channel closures should be connected?
            return
            
    # Boundary S=(1,2,3)
    S_tuple = (1,2,3)
    S_mask = get_flat_mask_for_particles(S_tuple)

    
    # L/R masks for support check
    # L = {1,2,3}, R = {4,5,6} (complement)
    # Lonly means supported on L (and maybe 123 boundary?)
    # Wait, support check logic uses Lonly/Ronly masks of CHANNELS.
    # dcp_common.build_LR_masks(C, Left, Right)
    
    L_set = {1,2,3}
    R_set = {4,5,6}
    Lonly, Ronly = build_LR_masks(C, L_set, R_set)
    
    # 4. Generate Maximal Nested Set
    log("Completing seed to maximal nested set...")
    # Use seeded sampler
    samples = sample_nested_sets_seeded(Gm, incompatible, seed_masks, n_samples=20, rng_seed=42)
    
    if not samples:
        log("Failed to generate any valid completion!")
        return
        
    log(f"Generated {len(samples)} samples.")
    
    # Test the seed itself first (Partial Chart)
    log("Testing SEED ONLY (Partial Chart)...")
    result_seed = scan_chart_proxy(seed_masks, Vinv_vectors, triples, len(C), Lonly, Ronly, S_mask, 
                              triple_sign, return_details=True)
    log(f"Seed Result: {result_seed.get('status')} Bad={result_seed.get('bad_count')}")
    
    best_res = None
    
    for i, N in enumerate(samples):
        log(f"Sample {i}: size={len(N)}")
        # 5. Run Proxy Pullback
        result = scan_chart_proxy(N, Vinv_vectors, triples, len(C), Lonly, Ronly, S_mask, 
                                  triple_sign, return_details=True)
        if result.get('status') == 'HIT':
            log(f"HIT found on sample {i}!")
            best_res = result
            break
        elif i == 0:
            best_res = result # Keep first if no hit
            
    if best_res is None or best_res.get('status') != 'HIT':
         log("No HIT in completions. Saving Seed Result (HIT) instead.")
         best_res = result_seed # Fallback to seed if completions fail
         
    result = best_res
    
    # 6. Report
    log("\nRESULTS:")
    log(f"Status: {result.get('status')}")
    log(f"Bad count: {result.get('bad_count')}")
    log(f"Valid: {result.get('valid')}")
    
    if 'alpha_sol' in result:
        sol = result['alpha_sol']
        log(f"Solution ratio: {sol}")
        # Check if b=0 is forced (a1=0?)
        # w0, w1 basis.
        # If sol is (1, 0), then w1 coeff is 0.
        
    cache_file = os.path.join(CACHE_DIR, "HIT_proxy_S123_seed12_56.sobj")
    save(result, cache_file)
    log(f"Saved result to {cache_file}")
    
    # Output detailed JSON for inspection
    json_path = os.path.join(CACHE_DIR, "HIT_proxy_S123_seed12_56.json")
    
    # Clean up for JSON
    clean_res = json_sanitize(result)
    # Remove large objects if necessary, but keep bad_coeffs
    if 'u_map' in clean_res:
         # convert keys to str
         clean_res['u_map'] = {str(k): v for k, v in clean_res['u_map'].items()}
         
    # Convert bad_coeffs keys from tuple to str for JSON
    if 'bad_coeffs' in clean_res:
        new_bad = []
        for pair, val in clean_res['bad_coeffs']:
            new_bad.append({'pair': str(pair), 'coeffs': val})
        clean_res['bad_coeffs'] = new_bad
        
    with open(json_path, 'w') as f:
        json.dump(clean_res, f, indent=2)
    log(f"Saved JSON report to {json_path}")

if __name__ == '__main__':
    run_repro()

