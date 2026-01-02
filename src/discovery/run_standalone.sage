#!/usr/bin/env sage
"""
Standalone Discovery Runner
============================

Runs all phases with minimal dependencies.
"""

from sage.all import *
import sys
import os
import json
import time
from datetime import datetime

# Setup paths
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, 'src'))

os.makedirs('discovery_results', exist_ok=True)

def log(msg):
    """Simple logging."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")
    with open('discovery_results/discovery.log', 'a') as f:
        f.write(f"[{timestamp}] {msg}\n")

log("="*70)
log("GRAVITY POSITIVE GEOMETRY DISCOVERY")
log("="*70)

# Import required modules
try:
    from src.posgeom.forest_polytope import get_forest_polynomial, get_forest_exponents
    from src.posgeom.physics_map import eval_edge_vars_from_spinors
    from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
    from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
    from src.chy_oracle.amplitude_spinor import ang_bracket
    log("Core modules imported successfully")
except Exception as e:
    log(f"ERROR importing core modules: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Try to import pushforward
try:
    from src.posgeom.saddle_pushforward import compute_pushforward_saddle, moment_map_and_jacobian
    HAS_PUSHFORWARD = True
    log("Pushforward module available")
except ImportError:
    HAS_PUSHFORWARD = False
    log("WARNING: Pushforward module not available")

# Run Phase 1
log("\n" + "="*70)
log("PHASE 1: Pushforward Diagnostic")
log("="*70)

phase1_result = {'success': False, 'error': 'not_run'}

if HAS_PUSHFORWARD:
    try:
        # Import phase 1 functions
        exec(open(os.path.join(project_root, 'src', 'discovery', 'phase1_pushforward.sage')).read())
        phase1_result = run_phase1_parallel(n_trials=1000, n_cores=4)
        log(f"Phase 1 complete: {phase1_result.get('verdict', 'unknown')}")
    except Exception as e:
        log(f"Phase 1 failed: {e}")
        import traceback
        traceback.print_exc()
        phase1_result = {'success': False, 'error': str(e)}
else:
    log("Phase 1 skipped: pushforward module not available")

# Run Phase 2
log("\n" + "="*70)
log("PHASE 2: BCFW Amplituhedron")
log("="*70)

phase2_result = {'success': False, 'error': 'not_run'}

try:
    exec(open(os.path.join(project_root, 'src', 'discovery', 'phase2_bcfw.sage')).read())
    phase2_result = run_phase2_parallel(n_trials=10000, n_cores=8)
    log(f"Phase 2 complete: Best formula = {phase2_result.get('best_formula', 'none')}")
except Exception as e:
    log(f"Phase 2 failed: {e}")
    import traceback
    traceback.print_exc()
    phase2_result = {'success': False, 'error': str(e)}

# Run Phase 3
log("\n" + "="*70)
log("PHASE 3: Positivity Search")
log("="*70)

phase3_result = {'success': False, 'error': 'not_run'}

try:
    exec(open(os.path.join(project_root, 'src', 'discovery', 'phase3_positivity.sage')).read())
    phase3_result = run_phase3_parallel(n_brute=10000, n_optimization_restarts=100, n_cores=4)
    log(f"Phase 3 complete: Positive region found = {phase3_result.get('positive_region_found', False)}")
except Exception as e:
    log(f"Phase 3 failed: {e}")
    import traceback
    traceback.print_exc()
    phase3_result = {'success': False, 'error': str(e)}

# Save all results
all_results = {
    'phase1': phase1_result,
    'phase2': phase2_result,
    'phase3': phase3_result,
    'timestamp': datetime.now().isoformat()
}

with open('discovery_results/all_phases_results.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=str)

# Generate verdict
log("\n" + "="*70)
log("GENERATING VERDICT")
log("="*70)

try:
    exec(open(os.path.join(project_root, 'src', 'discovery', 'verdict.sage')).read())
    verdict_result = main()
    log(f"\nVERDICT: {verdict_result.get('verdict', 'unknown')}")
    log(f"Confidence: {verdict_result.get('confidence', 'unknown')}")
except Exception as e:
    log(f"Verdict generation failed: {e}")
    import traceback
    traceback.print_exc()

log("\n" + "="*70)
log("DISCOVERY COMPLETE")
log("="*70)
log("Check discovery_results/VERDICT.md for final conclusion")


