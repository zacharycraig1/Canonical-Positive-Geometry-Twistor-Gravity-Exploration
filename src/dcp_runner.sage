# dcp_runner.sage
import sys, os
from sage.all import *

# Load modules
load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_invariants.sage")
load("src/dcp_search.sage")
load("src/dcp_intersect.sage")

def main():
    global INVARIANT_MODE  # Declare global at the start of the function
    # Allow overriding INVARIANT_MODE via global before running
    # but here we rely on dcp_common defaults or what we set here
    
    t_all = time.time()
    
    print("=" * 70)
    print("DCP GRAVITY - MODULAR RUNNER")
    print("=" * 70)
    log("Optimizations applied: PARI kernel, Scaled workers, Modular prefilter")
    log(f"Workers: {SCAN_WORKERS} (detected {_TOTAL_CORES} cores)")
    
    if FORCE_FRESH:
        clear_all_caches()
    
    # Initialize sign table
    sign_table = get_sign_table()
    
    # Build matroid
    log("\nBuilding M6 matroid...")
    C, M6 = build_M6_matroid()
    m = len(C)
    log(f"  {m} channels, rank {M6.rank()}")
    
    # Build OS3
    log("\nBuilding OS3 space...")
    triples, Vbasis = build_OS3_data(C, M6)
    log(f"  OS3 dim = {Vbasis.nrows()}")
    
    triples_array = np.array(triples, dtype=np.int32)
    
    # Invariants (initial computation)
    log("Computing invariants...")
    
    # We will let the loop handle invariant computation if multi-strategy is on
    Vinv = None 
    if not MULTI_STRATEGY_SEARCH:
        if INVARIANT_MODE == 'S6':
            Vinv = compute_S6_invariants(C, triples, Vbasis)
        elif INVARIANT_MODE in ('S3xS3','S3xS3Z2'):
            Vinv = compute_boundary_stabilizer_invariants(C, triples, Vbasis, mode=INVARIANT_MODE)
        else:
            raise ValueError(f"Unknown INVARIANT_MODE={INVARIANT_MODE}")
        log(f"  {INVARIANT_MODE} invariants: {len(Vinv)} vectors")
    
    # Load flats and incompatibility
    log("\nLoading connected flats...")
    Gm, full_mask = get_connected_flats(M6)
    Gset = set(int(x) for x in Gm)

    log("\nLoading incompatibility...")
    incompatible = get_incompatible_pairs(M6, Gm, Gset)
    
    if INTERSECTION_MODE:
        if MULTI_STRATEGY_SEARCH:
            log("\n" + "="*70)
            log("MULTI-STRATEGY SEARCH ENABLED")
            log(f"Will try invariant modes: {STRATEGY_INVARIANT_MODES}")
            log("="*70)
            
            best_result = None
            best_mode = None
            best_dim = None
            
            # We need to access global INVARIANT_MODE to set it for reporting
            # global INVARIANT_MODE # Moved to top
            original_mode_setting = INVARIANT_MODE
            
            for strategy_idx, mode in enumerate(STRATEGY_INVARIANT_MODES):
                log(f"\n{'='*70}")
                log(f"STRATEGY {strategy_idx+1}/{len(STRATEGY_INVARIANT_MODES)}: {mode}")
                log(f"{'='*70}")
                
                # Recompute invariants for this mode
                if mode == 'S6':
                    Vinv_strategy = compute_S6_invariants(C, triples, Vbasis)
                elif mode in ('S3xS3', 'S3xS3Z2'):
                    Vinv_strategy = compute_boundary_stabilizer_invariants(C, triples, Vbasis, mode=mode)
                else:
                    log(f"  Unknown mode {mode}, skipping")
                    continue
                
                log(f"  {mode} invariants: {len(Vinv_strategy)} vectors")
                
                # Run intersection with this strategy
                alpha0_strategy, table_strategy, Tmat_strategy, Vinv_final_strategy = intersection_across_boundaries(
                    C, triples, triples_array, m, Gm, Gset, full_mask, incompatible, 
                    Vinv_strategy, Vbasis, triple_sign
                )
                
                final_dim_strategy = len(Vinv_final_strategy) if Vinv_final_strategy is not None else None
                
                log(f"\n  Strategy {mode} result: dim={final_dim_strategy}")
                
                # Check if this is better
                if alpha0_strategy is not None:
                    if best_result is None:
                        best_result = (alpha0_strategy, table_strategy, Tmat_strategy, Vinv_final_strategy)
                        best_mode = mode
                        best_dim = final_dim_strategy
                    elif final_dim_strategy == 1 and best_dim != 1:
                        best_result = (alpha0_strategy, table_strategy, Tmat_strategy, Vinv_final_strategy)
                        best_mode = mode
                        best_dim = final_dim_strategy
                        log(f"  ✓ New best result: dim={final_dim_strategy} (mode={mode})")
                        
                        if STRATEGY_STOP_ON_DIM1:
                            log(f"\n  ✓ Found dim=1 candidate with mode {mode}! Stopping search.")
                            break
                    elif final_dim_strategy < best_dim:
                        best_result = (alpha0_strategy, table_strategy, Tmat_strategy, Vinv_final_strategy)
                        best_mode = mode
                        best_dim = final_dim_strategy
                        log(f"  ✓ New best result: dim={final_dim_strategy} (mode={mode})")
                else:
                    log(f"  ✗ No candidate found with mode {mode}")
            
            if best_result is None:
                log("\n[STRATEGY] No candidate found with any invariant mode.")
                return
            
            log(f"\n{'='*70}")
            log(f"BEST RESULT: mode={best_mode}, dim={best_dim}")
            log(f"{'='*70}")
            
            alpha0, table, Tmat, Vinv_final = best_result
            INVARIANT_MODE = best_mode
        else:
            # Single strategy
            alpha0, table, Tmat, Vinv_final = intersection_across_boundaries(
                C, triples, triples_array, m, Gm, Gset, full_mask, incompatible, Vinv, Vbasis, triple_sign
            )
        
        # Report & Verify
        if alpha0 is None:
            log("\n[INTERSECT] Empty intersection (no surviving candidate).")
            return

        final_dim = len(Vinv_final)
        
        # Artifact saving
        run_id = time.strftime("%Y-%m-%d_%H%M%S") + f"_seed{int(BASE_SEED)}_mode{str(INVARIANT_MODE)}"
        run_dir = os.path.join(os.getcwd(), ARTIFACT_DIR, run_id)
        if SAVE_ARTIFACTS:
            os.makedirs(run_dir, exist_ok=True)
            
            # Save artifacts
            save(Tmat, os.path.join(run_dir, "candidate_space_basis_in_invariant_coords.sobj"))
            save(Vinv_final, os.path.join(run_dir, "candidate_space_basis_in_OS_coords.sobj"))
            save(alpha0, os.path.join(run_dir, "alpha0_first_basis_vector.sobj"))
            
            report = {
                "mode": "intersection",
                "invariant_mode": str(INVARIANT_MODE),
                "final_dim": int(final_dim),
                "intersection_log": table,
            }
            with open(os.path.join(run_dir, "hit_report_intersection.json"), "w") as f:
                json.dump(json_sanitize(report), f, indent=2)
                
            log(f"[INTERSECT] Wrote report to {run_dir}")

        # Verification if dim=1
        if final_dim == 1:
            log("[CANDIDATE] Candidate direction found.")
            verify_path = os.path.join(run_dir, "dim1_candidate_verification.json") if SAVE_ARTIFACTS else None
            verification = verify_dim1_candidate(
                alpha0, Vinv_final, Tmat, C, triples, Vbasis,
                check_hodge=True,
                save_path=verify_path
            )
            if verification.get('status') == 'verified':
                log("\n[SUCCESS] Dimension-1 candidate verified!")
        else:
            log(f"[CANDIDATE] Not unique yet (dim = {final_dim}).")

    log(f"Total time: {time.time() - t_all:.1f}s")

if __name__ == "__main__":
    main()

