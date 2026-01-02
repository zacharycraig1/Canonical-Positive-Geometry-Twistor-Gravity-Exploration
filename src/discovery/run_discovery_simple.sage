#!/usr/bin/env sage
"""
Simplified Discovery Runner
============================

Runs phases sequentially to avoid import issues.
"""

from sage.all import *
import sys
import os

# Add project root
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, 'src'))

# Load modules using Sage's load mechanism
load(os.path.join(project_root, 'src', 'discovery', 'utils.sage'))
load(os.path.join(project_root, 'src', 'discovery', 'phase1_pushforward.sage'))
load(os.path.join(project_root, 'src', 'discovery', 'phase2_bcfw.sage'))
load(os.path.join(project_root, 'src', 'discovery', 'phase3_positivity.sage'))
load(os.path.join(project_root, 'src', 'discovery', 'verdict.sage'))

# Now run
print("="*70)
print("GRAVITY POSITIVE GEOMETRY DISCOVERY")
print("="*70)
print()

# Run phases sequentially
results = {}

print("Starting Phase 1: Pushforward Diagnostic...")
try:
    results['phase1'] = run_phase1_parallel(n_trials=1000, n_cores=4)
    print(f"Phase 1 complete: {results['phase1'].get('verdict', 'unknown')}")
except Exception as e:
    print(f"Phase 1 failed: {e}")
    results['phase1'] = {'success': False, 'error': str(e)}

print("\nStarting Phase 2: BCFW Amplituhedron...")
try:
    results['phase2'] = run_phase2_parallel(n_trials=10000, n_cores=8)
    print(f"Phase 2 complete: Best formula = {results['phase2'].get('best_formula', 'none')}")
except Exception as e:
    print(f"Phase 2 failed: {e}")
    results['phase2'] = {'success': False, 'error': str(e)}

print("\nStarting Phase 3: Positivity Search...")
try:
    results['phase3'] = run_phase3_parallel(n_brute=10000, n_optimization_restarts=100, n_cores=4)
    print(f"Phase 3 complete: Positive region found = {results['phase3'].get('positive_region_found', False)}")
except Exception as e:
    print(f"Phase 3 failed: {e}")
    results['phase3'] = {'success': False, 'error': str(e)}

# Save results
os.makedirs('discovery_results', exist_ok=True)
import json
with open('discovery_results/all_phases_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\nGenerating verdict...")
try:
    verdict_result = main()  # From verdict.sage
    print("\n" + "="*70)
    print("DISCOVERY COMPLETE")
    print("="*70)
    print(f"Verdict: {verdict_result.get('verdict', 'unknown')}")
    print(f"Confidence: {verdict_result.get('confidence', 'unknown')}")
    print("="*70)
except Exception as e:
    print(f"Verdict generation failed: {e}")
    import traceback
    traceback.print_exc()


