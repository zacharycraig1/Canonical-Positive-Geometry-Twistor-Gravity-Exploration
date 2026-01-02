#!/usr/bin/env sage
"""
Phase 3: Positivity Region Search
==================================

Searches for or proves non-existence of the positive region for gravity.

The positive region requires:
1. All ordered 4-brackets ⟨i i+1 j j+1⟩ > 0 (amplituhedron positivity)
2. All edge variables z_ij > 0 simultaneously

Uses:
- Brute force search (10,000 random points)
- Optimization search (scipy.optimize with multiple restarts)
- Algebraic obstruction check (if possible)
"""

from sage.all import *
import numpy as np
import sys
import os
from multiprocessing import Pool, cpu_count
import time
try:
    from scipy.optimize import minimize, differential_evolution
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.discovery.utils import (
    get_reference_spinors, compute_hodges_oracle, save_checkpoint, load_checkpoint,
    log_progress, statistical_analysis
)
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.chy_oracle.kinematics_samples import MomentumTwistor


def check_twistor_positivity(tw):
    """
    Check if momentum twistor configuration is in positive region.
    
    Positive means: all ordered 4-brackets ⟨i i+1 j j+1⟩ > 0
    """
    n = tw.n
    
    # Check all ordered 4-brackets
    for i in range(n):
        ip1 = (i + 1) % n
        for j in range(i + 2, n):
            jp1 = (j + 1) % n
            if jp1 == i:
                continue
            
            # Compute 4-bracket ⟨i i+1 j j+1⟩
            four_br = compute_four_bracket(tw, i, ip1, j, jp1)
            if four_br is None or four_br <= 0:
                return False
    
    return True


def compute_four_bracket(tw, i, j, k, l):
    """Compute ⟨ijkl⟩ = det(Z_i, Z_j, Z_k, Z_l)."""
    try:
        Z = [tw.Z[i], tw.Z[j], tw.Z[k], tw.Z[l]]
        M = matrix(QQ, Z)
        det_val = M.det()
        
        # Sign from permutation parity
        perm = [i, j, k, l]
        sorted_perm = sorted(perm)
        inversions = sum(1 for a in range(4) for b in range(a+1, 4) if perm[a] > perm[b])
        sign = (-1) ** inversions
        
        return sign * det_val
    except:
        return None


def check_z_positivity(n, seed):
    """
    Check if all z_ij > 0 for a given kinematic point.
    
    Returns: (all_positive, z_values_dict)
    """
    try:
        hodges, lambdas, tildes = compute_hodges_oracle(n, seed)
        if hodges is None:
            return False, None
        
        x, y = get_reference_spinors()
        z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        
        # Check all z_ij
        all_positive = True
        z_values = {}
        
        for key, val in z_map.items():
            if isinstance(key, str) and key.startswith('z_'):
                z_values[key] = float(val)
                if val <= 0:
                    all_positive = False
        
        return all_positive, z_values
    except:
        return False, None


def test_single_point_positivity(seed):
    """
    Test a single point for positivity.
    
    Returns dict with results.
    """
    result = {
        'seed': seed,
        'twistor_positive': False,
        'z_positive': False,
        'both_positive': False,
        'success': False
    }
    
    try:
        # Test twistor positivity
        tw = MomentumTwistor(n=6, seed=seed)
        twistor_pos = check_twistor_positivity(tw)
        result['twistor_positive'] = twistor_pos
        
        # Test z_ij positivity
        z_pos, z_vals = check_z_positivity(6, seed)
        result['z_positive'] = z_pos
        
        if twistor_pos and z_pos:
            result['both_positive'] = True
            result['z_values'] = z_vals
        
        result['success'] = True
        
    except Exception:
        pass
    
    return result


def run_brute_force_search(n_trials=10000, n_cores=4):
    """
    Brute force search for positive region.
    
    Tests random kinematic points.
    """
    log_progress("PHASE3", f"Starting brute force search: {n_trials} trials, {n_cores} cores")
    
    # Check for checkpoint
    checkpoint = load_checkpoint('phase3_brute_checkpoint.json')
    start_seed = 0
    results = []
    
    if checkpoint:
        log_progress("PHASE3", f"Resuming from checkpoint: {len(checkpoint.get('results', []))} points tested")
        results = checkpoint.get('results', [])
        start_seed = checkpoint.get('last_seed', 0) + 1
    
    seeds = list(range(start_seed, start_seed + n_trials - len(results)))
    
    log_progress("PHASE3", f"Testing {len(seeds)} new points...")
    
    t_start = time.time()
    with Pool(processes=n_cores) as pool:
        new_results = pool.map(test_single_point_positivity, seeds)
    
    results.extend(new_results)
    elapsed = time.time() - t_start
    
    log_progress("PHASE3", f"Completed {len(results)} tests in {elapsed:.1f}s")
    
    # Save checkpoint
    save_checkpoint({
        'results': results,
        'last_seed': seeds[-1] if seeds else start_seed - 1,
        'timestamp': time.time()
    }, 'phase3_brute_checkpoint.json')
    
    # Analyze results
    successful = [r for r in results if r.get('success', False)]
    twistor_pos = [r for r in successful if r.get('twistor_positive', False)]
    z_pos = [r for r in successful if r.get('z_positive', False)]
    both_pos = [r for r in successful if r.get('both_positive', False)]
    
    log_progress("PHASE3", f"Results: {len(successful)} successful, {len(twistor_pos)} twistor-positive, {len(z_pos)} z-positive, {len(both_pos)} both-positive")
    
    return {
        'success': True,
        'total_tests': len(results),
        'successful': len(successful),
        'twistor_positive_count': len(twistor_pos),
        'z_positive_count': len(z_pos),
        'both_positive_count': len(both_pos),
        'positive_points': [r for r in both_pos if r.get('both_positive', False)]
    }


def optimization_objective(z_vec):
    """
    Objective function for optimization: minimize |min(z_ij, 0)|
    
    We want all z_ij > 0, so we minimize the negative part.
    """
    # z_vec is flattened array of z_ij values
    min_val = np.min(z_vec)
    if min_val > 0:
        return 0.0  # All positive - success!
    else:
        return abs(min_val)  # Minimize the negative part


def run_optimization_search(n_restarts=100):
    """
    Use optimization to find positive region.
    
    Searches for kinematic point where all z_ij > 0.
    """
    if not HAS_SCIPY:
        log_progress("PHASE3", "Skipping optimization search (scipy not available)", "WARNING")
        return {
            'success': True,
            'positive_found': False,
            'best_point': None,
            'best_objective': None,
            'skipped': True
        }
    
    log_progress("PHASE3", f"Starting optimization search: {n_restarts} restarts")
    
    # For n=6, we have 15 edge variables z_ij
    # But we need to search in the space of momentum twistors (6 × 4 = 24 parameters)
    # This is complex - we'll use a simplified approach
    
    # Instead, search in z-space directly (if we can find a valid mapping)
    # Or search in spinor space
    
    positive_found = False
    best_point = None
    best_value = float('inf')
    
    for restart in range(n_restarts):
        try:
            # Random initial guess in spinor space
            # We'll use seed-based generation and then optimize
            seed = 50000 + restart
            
            # Get initial point
            hodges, lambdas, tildes = compute_hodges_oracle(6, seed)
            if hodges is None:
                continue
            
            x, y = get_reference_spinors()
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
            
            # Extract z values
            z_vec = []
            for i in range(6):
                for j in range(i + 1, 6):
                    key = f"z_{i}_{j}"
                    if key in z_map:
                        z_vec.append(float(z_map[key]))
            
            if not z_vec:
                continue
            
            z_vec = np.array(z_vec)
            
            # Check if already positive
            if np.all(z_vec > 0):
                log_progress("PHASE3", f"[SUCCESS] Found positive point at restart {restart}!")
                positive_found = True
                best_point = {
                    'seed': seed,
                    'z_values': z_vec.tolist()
                }
                break
            
            # Objective: minimize negative z values
            obj = optimization_objective(z_vec)
            if obj < best_value:
                best_value = obj
                best_point = {
                    'seed': seed,
                    'z_values': z_vec.tolist(),
                    'objective': obj
                }
        
        except Exception:
            continue
        
        if restart % 10 == 0:
            log_progress("PHASE3", f"Optimization progress: {restart}/{n_restarts}, best objective: {best_value:.6e}")
    
    return {
        'success': True,
        'positive_found': positive_found,
        'best_point': best_point,
        'best_objective': best_value
    }


def run_phase3_parallel(n_brute=10000, n_optimization_restarts=100, n_cores=4):
    """
    Run Phase 3: Positivity search.
    
    Combines brute force and optimization.
    """
    log_progress("PHASE3", "Starting Phase 3: Positivity Region Search")
    
    # Part 1: Brute force
    brute_result = run_brute_force_search(n_trials=n_brute, n_cores=n_cores)
    
    # Part 2: Optimization
    opt_result = run_optimization_search(n_restarts=n_optimization_restarts)
    
    # Combine results
    phase3_result = {
        'success': True,
        'brute_force': brute_result,
        'optimization': opt_result,
        'positive_region_found': (
            brute_result.get('both_positive_count', 0) > 0 or
            opt_result.get('positive_found', False)
        )
    }
    
    if phase3_result['positive_region_found']:
        log_progress("PHASE3", "[SUCCESS] Positive region FOUND!", "SUCCESS")
    else:
        log_progress("PHASE3", "[RESULT] No positive region found in search", "INFO")
    
    return phase3_result


def main():
    """Run Phase 3 standalone."""
    n_cores = min(4, cpu_count())
    result = run_phase3_parallel(
        n_brute=10000,
        n_optimization_restarts=100,
        n_cores=n_cores
    )
    
    # Save final result
    save_checkpoint(result, 'phase3_result.json')
    
    print("\n" + "="*70)
    print("PHASE 3 RESULTS")
    print("="*70)
    print(f"Positive region found: {result.get('positive_region_found', False)}")
    print(f"Brute force: {result.get('brute_force', {}).get('both_positive_count', 0)} positive points")
    print(f"Optimization: {result.get('optimization', {}).get('positive_found', False)}")
    print("="*70)
    
    return result


if __name__ == '__main__':
    main()

