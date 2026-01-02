#!/usr/bin/env sage
"""
Phase 1: Pushforward Diagnostic
================================

Systematic analysis of why saddle pushforward works for n=4,5 but fails for n=6.

Tests 1000 kinematic points to determine:
1. Is pushforward/Hodges ratio constant? (→ normalization fix possible)
2. Is ratio varying? (→ fundamental problem with map)
3. What is the Jacobian structure?
"""

from sage.all import *
import numpy as np
import sys
import os
from multiprocessing import Pool, cpu_count
import time

# Add project root
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, 'src'))

# Import utils - handle both .sage and .py
try:
    from src.discovery.utils import (
        get_reference_spinors, compute_hodges_oracle, compute_forest_polynomial_value,
        is_singular_kinematics, save_checkpoint, load_checkpoint, log_progress,
        statistical_analysis, is_constant_ratio
    )
except ImportError:
    # Fallback: define minimal versions
    def get_reference_spinors():
        return vector(QQ, [1, 0]), vector(QQ, [0, 1])
    
    def compute_hodges_oracle(n, seed, roots=(0, 1, 2)):
        from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
        from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
        try:
            lambdas, tildes = sample_spinors_from_twistor(n=n, seed=seed)
            x, y = get_reference_spinors()
            M, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=roots)
            return M, lambdas, tildes
        except:
            return None, None, None
    
    def log_progress(phase, message, level="INFO"):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] [{phase}] [{level}] {message}")
    
    def save_checkpoint(data, filename):
        os.makedirs('discovery_results', exist_ok=True)
        import json
        with open(os.path.join('discovery_results', filename), 'w') as f:
            json.dump(data, f, indent=2, default=str)
    
    def load_checkpoint(filename):
        import json
        filepath = os.path.join('discovery_results', filename)
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                return json.load(f)
        return None
    
    def statistical_analysis(values, name="values"):
        if not values:
            return None
        values_float = [float(v) for v in values if v is not None]
        if not values_float:
            return None
        mean_val = sum(values_float) / len(values_float)
        variance = sum((x - mean_val)**2 for x in values_float) / len(values_float)
        std_val = variance ** 0.5
        return {
            'mean': mean_val,
            'std': std_val,
            'min': min(values_float),
            'max': max(values_float),
            'cv': std_val / abs(mean_val) if mean_val != 0 else float('inf'),
            'count': len(values_float)
        }
    
    def is_constant_ratio(ratios, threshold=0.01):
        stats = statistical_analysis(ratios)
        if stats is None:
            return False, None
        is_const = stats['cv'] < threshold
        return is_const, stats
    
    from datetime import datetime
from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.physics_map import eval_edge_vars_from_spinors

# Try to import pushforward
try:
    from src.posgeom.saddle_pushforward import compute_pushforward_saddle, moment_map_and_jacobian
    HAS_PUSHFORWARD = True
except ImportError:
    HAS_PUSHFORWARD = False
    log_progress("PHASE1", "Warning: saddle_pushforward not available", "WARNING")


def test_single_point(seed, n=6, roots=[0, 1, 2]):
    """
    Test a single kinematic point.
    
    Returns dict with:
    - seed: seed used
    - hodges: Hodges value
    - pushforward: pushforward value
    - ratio: pushforward/hodges
    - forest_val: forest polynomial value
    - success: bool
    """
    result = {
        'seed': seed,
        'hodges': None,
        'pushforward': None,
        'ratio': None,
        'forest_val': None,
        'success': False
    }
    
    try:
        # Get Hodges oracle
        hodges, lambdas, tildes = compute_hodges_oracle(n, seed, roots=tuple(roots))
        if hodges is None or hodges == 0:
            return result
        
        result['hodges'] = float(hodges)
        
        # Check singularity
        if is_singular_kinematics(lambdas, tildes):
            return result
        
        # Get forest polynomial exponents
        exponents, edge_order = get_forest_exponents(n, roots)
        poly_exponents = np.array(exponents)
        poly_coeffs = np.ones(len(exponents))
        
        # Compute z_ij values
        x, y = get_reference_spinors()
        z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        
        # Get z values in canonical edge order
        z_vec = []
        for edge in edge_order:
            i, j = edge
            key = f"z_{i}_{j}"
            if key not in z_map:
                return result
            z_val = z_map[key]
            if z_val <= 0:
                return result  # Outside positive region
            z_vec.append(float(z_val))
        
        z_vec = np.array(z_vec)
        log_z = np.log(z_vec)
        
        # Compute moment map X = grad_logz log(F)
        X, J = moment_map_and_jacobian(log_z, poly_coeffs, poly_exponents)
        det_J = np.linalg.det(J)
        
        if abs(det_J) < 1e-15:
            return result
        
        pushforward_val = 1.0 / det_J
        result['pushforward'] = float(pushforward_val)
        
        # Compute ratio
        if hodges != 0:
            ratio = pushforward_val / float(hodges)
            result['ratio'] = float(ratio)
        
        # Forest polynomial value
        forest_val = compute_forest_polynomial_value(n, roots, lambdas, tildes, x, y)
        if forest_val is not None:
            result['forest_val'] = float(forest_val)
        
        result['success'] = True
        
    except Exception as e:
        pass
    
    return result


def run_phase1_parallel(n_trials=1000, n_cores=4):
    """
    Run Phase 1 tests in parallel.
    
    Returns dict with all results and statistics.
    """
    log_progress("PHASE1", f"Starting Phase 1: Pushforward Diagnostic ({n_trials} trials, {n_cores} cores)")
    
    if not HAS_PUSHFORWARD:
        log_progress("PHASE1", "ERROR: Pushforward module not available", "ERROR")
        return {
            'success': False,
            'error': 'pushforward_module_not_available'
        }
    
    # Check for checkpoint
    checkpoint = load_checkpoint('phase1_checkpoint.json')
    start_seed = 0
    results = []
    
    if checkpoint:
        log_progress("PHASE1", f"Resuming from checkpoint: {len(checkpoint.get('results', []))} points already tested")
        results = checkpoint.get('results', [])
        start_seed = checkpoint.get('last_seed', 0) + 1
    
    # Generate seed range
    seeds = list(range(start_seed, start_seed + n_trials - len(results)))
    
    log_progress("PHASE1", f"Testing {len(seeds)} new points...")
    
    # Run in parallel
    t_start = time.time()
    with Pool(processes=n_cores) as pool:
        new_results = pool.starmap(test_single_point, [(seed,) for seed in seeds])
    
    results.extend(new_results)
    elapsed = time.time() - t_start
    
    log_progress("PHASE1", f"Completed {len(results)} tests in {elapsed:.1f}s")
    
    # Save checkpoint
    save_checkpoint({
        'results': results,
        'last_seed': seeds[-1] if seeds else start_seed - 1,
        'timestamp': time.time()
    }, 'phase1_checkpoint.json')
    
    # Analyze results
    successful = [r for r in results if r['success']]
    log_progress("PHASE1", f"Successful tests: {len(successful)}/{len(results)}")
    
    if len(successful) < 10:
        log_progress("PHASE1", "WARNING: Too few successful tests for statistical analysis", "WARNING")
        return {
            'success': False,
            'error': 'insufficient_data',
            'successful_count': len(successful)
        }
    
    # Extract ratios
    ratios = [r['ratio'] for r in successful if r['ratio'] is not None]
    
    if not ratios:
        log_progress("PHASE1", "ERROR: No valid ratios computed", "ERROR")
        return {
            'success': False,
            'error': 'no_ratios'
        }
    
    # Statistical analysis
    stats = statistical_analysis(ratios, "pushforward_ratios")
    is_const, const_stats = is_constant_ratio(ratios, threshold=0.01)
    
    log_progress("PHASE1", f"Ratio statistics: mean={stats['mean']:.6e}, std={stats['std']:.6e}, cv={stats['cv']:.6e}")
    
    if is_const:
        log_progress("PHASE1", f"[SUCCESS] Ratio is CONSTANT! Normalization factor: {1/stats['mean']:.6e}", "SUCCESS")
        verdict = "constant_ratio"
        normalization = 1.0 / stats['mean']
    else:
        log_progress("PHASE1", "[RESULT] Ratio is NOT constant - map may be fundamentally wrong", "INFO")
        verdict = "varying_ratio"
        normalization = None
    
    # Jacobian analysis (sample)
    jacobian_dets = []
    for r in successful[:100]:  # Sample first 100
        try:
            seed = r['seed']
            hodges, lambdas, tildes = compute_hodges_oracle(6, seed)
            if hodges is None:
                continue
            
            x, y = get_reference_spinors()
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
            exponents, edge_order = get_forest_exponents(6, [0, 1, 2])
            
            z_vec = []
            for edge in edge_order:
                i, j = edge
                z_vec.append(float(z_map[f"z_{i}_{j}"]))
            
            z_vec = np.array(z_vec)
            if np.any(z_vec <= 0):
                continue
            
            log_z = np.log(z_vec)
            poly_exponents = np.array(exponents)
            poly_coeffs = np.ones(len(exponents))
            
            X, J = moment_map_and_jacobian(log_z, poly_coeffs, poly_exponents)
            det_J = np.linalg.det(J)
            jacobian_dets.append(float(det_J))
        except:
            continue
    
    jacobian_stats = statistical_analysis(jacobian_dets, "jacobian_determinants")
    
    # Final result
    phase1_result = {
        'success': True,
        'verdict': verdict,
        'normalization': normalization,
        'ratio_stats': stats,
        'is_constant': is_const,
        'successful_count': len(successful),
        'total_tests': len(results),
        'jacobian_stats': jacobian_stats,
        'ratios_sample': ratios[:20]  # First 20 for inspection
    }
    
    log_progress("PHASE1", f"Phase 1 complete: {verdict}")
    
    return phase1_result


def main():
    """Run Phase 1 standalone."""
    n_cores = min(4, cpu_count())
    result = run_phase1_parallel(n_trials=1000, n_cores=n_cores)
    
    # Save final result
    save_checkpoint(result, 'phase1_result.json')
    
    print("\n" + "="*70)
    print("PHASE 1 RESULTS")
    print("="*70)
    print(f"Verdict: {result.get('verdict', 'unknown')}")
    if result.get('normalization'):
        print(f"Normalization factor: {result['normalization']:.6e}")
    print(f"Successful tests: {result.get('successful_count', 0)}")
    print("="*70)
    
    return result


if __name__ == '__main__':
    main()

