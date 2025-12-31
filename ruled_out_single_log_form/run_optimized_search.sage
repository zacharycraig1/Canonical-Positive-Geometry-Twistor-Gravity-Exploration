#!/usr/bin/env sage
# Optimized iterative search with automatic parameter adjustment
# Run with: sage run_optimized_search.sage

load('54.sage')
import json
import os
import sys

def load_suggestions():
    """Load suggestions from previous run if available."""
    runs_dir = "runs"
    if not os.path.exists(runs_dir):
        return None
    
    run_dirs = sorted([d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))], reverse=True)
    for run_dir in run_dirs[:3]:  # Check last 3 runs
        sugg_file = os.path.join(runs_dir, run_dir, "next_config_suggestions.json")
        if os.path.exists(sugg_file):
            try:
                with open(sugg_file, 'r') as f:
                    return json.load(f)
            except:
                continue
    return None

def apply_suggestions(suggestions):
    """Apply configuration suggestions."""
    global INVARIANT_MODE, TOTAL_TRIALS, BASE_SEED, INTERSECT_BOUNDARY_MODE, INTERSECT_BOUNDARIES
    
    if not suggestions:
        return
    
    if 'INVARIANT_MODE' in suggestions:
        INVARIANT_MODE = suggestions['INVARIANT_MODE']
        log(f"Applied: INVARIANT_MODE = {INVARIANT_MODE}")
    
    if 'TOTAL_TRIALS' in suggestions:
        TOTAL_TRIALS = suggestions['TOTAL_TRIALS']
        log(f"Applied: TOTAL_TRIALS = {TOTAL_TRIALS}")
    
    if 'BASE_SEED' in suggestions:
        BASE_SEED = suggestions['BASE_SEED']
        log(f"Applied: BASE_SEED = {BASE_SEED}")
    
    if 'INTERSECT_BOUNDARY_MODE' in suggestions:
        INTERSECT_BOUNDARY_MODE = suggestions['INTERSECT_BOUNDARY_MODE']
        log(f"Applied: INTERSECT_BOUNDARY_MODE = {INTERSECT_BOUNDARY_MODE}")
    
    if 'INTERSECT_BOUNDARIES' in suggestions:
        INTERSECT_BOUNDARIES = [tuple(b) for b in suggestions['INTERSECT_BOUNDARIES']]
        log(f"Applied: INTERSECT_BOUNDARIES = {INTERSECT_BOUNDARIES}")

def check_for_success(run_dir):
    """Check if we found dim=1 candidate."""
    report_file = os.path.join(run_dir, "hit_report_intersection.json")
    if not os.path.exists(report_file):
        return False
    
    try:
        with open(report_file, 'r') as f:
            report = json.load(f)
        return report.get('final_dim') == 1
    except:
        return False

def main_optimized():
    """Run optimized iterative search."""
    max_iterations = 10
    iteration = 1
    
    log("\n" + "="*70)
    log("OPTIMIZED ITERATIVE SEARCH")
    log("="*70)
    
    # Load suggestions from previous run
    suggestions = load_suggestions()
    if suggestions:
        log(f"\nLoaded suggestions from previous run:")
        for change in suggestions.get('changes', []):
            log(f"  - {change}")
        apply_suggestions(suggestions)
    
    while iteration <= max_iterations:
        log(f"\n{'='*70}")
        log(f"ITERATION {iteration}/{max_iterations}")
        log(f"{'='*70}")
        
        try:
            main()
        except KeyboardInterrupt:
            log("\nInterrupted by user")
            break
        except Exception as e:
            log(f"\nError in iteration {iteration}: {e}")
            import traceback
            traceback.print_exc()
            iteration += 1
            continue
        
        # Check for success
        runs_dir = "runs"
        if os.path.exists(runs_dir):
            run_dirs = sorted([d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))], reverse=True)
            if run_dirs:
                latest_run = os.path.join(runs_dir, run_dirs[0])
                if check_for_success(latest_run):
                    log("\n" + "!"*70)
                    log("SUCCESS! Found dim=1 candidate!")
                    log(f"Results in: {latest_run}")
                    log("!"*70)
                    return
        
        # Load new suggestions for next iteration
        suggestions = load_suggestions()
        if suggestions and suggestions.get('action') != 'SUCCESS - verify candidate':
            apply_suggestions(suggestions)
        
        iteration += 1
    
    log(f"\n{'='*70}")
    log("Reached max iterations")
    log(f"{'='*70}")

if __name__ == '__main__':
    main_optimized()




