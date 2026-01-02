#!/usr/bin/env sage
"""
Phase 2: BCFW Amplituhedron Verification
=========================================

Tests if BCFW cell sum equals Hodges determinant for gravity MHV.

Implements multiple BCFW formula variants:
1. Standard BCFW recursion
2. Gravity-doubled Parke-Taylor structure
3. KLT-inspired bilinear form

Tests 10,000 kinematic points to find the correct formula.
"""

from sage.all import *
import numpy as np
import sys
import os
from multiprocessing import Pool, cpu_count
import time
from itertools import combinations

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.discovery.utils import (
    get_reference_spinors, compute_hodges_oracle, save_checkpoint, load_checkpoint,
    log_progress, statistical_analysis, is_constant_ratio
)
from src.chy_oracle.kinematics_samples import MomentumTwistor


class BCFWFormula:
    """Base class for BCFW formula variants."""
    
    def __init__(self, name):
        self.name = name
    
    def compute(self, tw):
        """Compute amplitude from momentum twistor. Returns None if fails."""
        raise NotImplementedError


class StandardBCFW(BCFWFormula):
    """Standard BCFW recursion for gravity MHV."""
    
    def __init__(self):
        super().__init__("StandardBCFW")
    
    def compute(self, tw):
        """Standard BCFW: sum over channels with Parke-Taylor squared."""
        n = tw.n
        if n != 6:
            return None
        
        # Parke-Taylor factor (squared for gravity)
        pt = self._parke_taylor(tw)
        if pt is None:
            return None
        
        pt_squared = pt * pt
        
        # Sum over 2-particle channels
        total = QQ(0)
        for i in range(n):
            ip1 = (i + 1) % n
            s = self._mandelstam(tw, [i, ip1])
            if s is None or s == 0:
                continue
            total += pt_squared / s
        
        return total
    
    def _parke_taylor(self, tw):
        """Compute Parke-Taylor: 1 / (⟨12⟩⟨23⟩...⟨n1⟩)"""
        n = tw.n
        result = QQ(1)
        for i in range(n):
            ip1 = (i + 1) % n
            ang = tw.get_angle(i, ip1)
            if ang == 0:
                return None
            result /= ang
        return result
    
    def _mandelstam(self, tw, indices):
        """Compute Mandelstam invariant s_{indices}."""
        if len(indices) == 2:
            i, j = indices
            ang = tw.get_angle(i, j)
            sq = self._get_square(tw, i, j)
            if sq is None:
                return None
            return ang * sq
        else:
            # Sum over pairs
            total = QQ(0)
            for a, i in enumerate(indices):
                for b in range(a + 1, len(indices)):
                    j = indices[b]
                    s_ij = self._mandelstam(tw, [i, j])
                    if s_ij is None:
                        return None
                    total += s_ij
            return total
    
    def _get_square(self, tw, i, j):
        """Compute [ij] from momentum twistors using tilde_lambda reconstruction."""
        try:
            tilde_i = tw.get_tilde_lambda(i)
            tilde_j = tw.get_tilde_lambda(j)
            if tilde_i is None or tilde_j is None:
                return None
            # [ij] = det(tilde_i, tilde_j)
            return tilde_i[0] * tilde_j[1] - tilde_i[1] * tilde_j[0]
        except Exception:
            return None


class GravityDoubledPT(BCFWFormula):
    """Gravity as squared Parke-Taylor with additional factors."""
    
    def __init__(self):
        super().__init__("GravityDoubledPT")
    
    def compute(self, tw):
        """Gravity = (PT)^2 × correction_factor."""
        n = tw.n
        if n != 6:
            return None
        
        # Parke-Taylor
        pt = self._parke_taylor(tw)
        if pt is None:
            return None
        
        # Try different correction factors
        # Factor 1: Product of all angle brackets
        angle_product = QQ(1)
        for i in range(n):
            ip1 = (i + 1) % n
            ang = tw.get_angle(i, ip1)
            if ang == 0:
                return None
            angle_product *= ang
        
        # Gravity = (PT)^2 / (angle_product)^2
        result = (pt * pt) / (angle_product * angle_product)
        
        return result
    
    def _parke_taylor(self, tw):
        """Parke-Taylor factor."""
        n = tw.n
        result = QQ(1)
        for i in range(n):
            ip1 = (i + 1) % n
            ang = tw.get_angle(i, ip1)
            if ang == 0:
                return None
            result /= ang
        return result


class KLTBilinear(BCFWFormula):
    """KLT-inspired bilinear form: sum over permutations."""
    
    def __init__(self):
        super().__init__("KLTBilinear")
    
    def compute(self, tw):
        """KLT structure: M = sum S[α|β] A(α) A(β)."""
        n = tw.n
        if n != 6:
            return None
        
        # For MHV, this simplifies
        # Try: M = sum_{perms} PT(perm)^2 / (kinematic factors)
        
        # Simplified: sum over cyclic permutations
        total = QQ(0)
        
        for shift in range(n):
            # Cyclic permutation
            perm = [(i + shift) % n for i in range(n)]
            
            # Parke-Taylor for this permutation
            pt = QQ(1)
            for i in range(n):
                a = perm[i]
                b = perm[(i + 1) % n]
                ang = tw.get_angle(a, b)
                if ang == 0:
                    break
                pt /= ang
            else:
                # No break - all angles non-zero
                total += pt * pt
        
        return total if total != 0 else None


def test_single_point_bcfw(seed, formulas):
    """
    Test all BCFW formulas on a single kinematic point.
    
    Returns dict with results for each formula.
    """
    result = {
        'seed': seed,
        'hodges': None,
        'formulas': {}
    }
    
    try:
        # Get Hodges oracle
        hodges, lambdas, tildes = compute_hodges_oracle(6, seed)
        if hodges is None or hodges == 0:
            return result
        
        result['hodges'] = float(hodges)
        
        # Create momentum twistor from same seed
        tw = MomentumTwistor(n=6, seed=seed)
        
        # Test each formula
        for formula in formulas:
            try:
                bcfw_val = formula.compute(tw)
                if bcfw_val is None:
                    result['formulas'][formula.name] = None
                    continue
                
                bcfw_float = float(bcfw_val)
                if hodges != 0:
                    ratio = bcfw_float / float(hodges)
                    result['formulas'][formula.name] = {
                        'value': bcfw_float,
                        'ratio': float(ratio)
                    }
                else:
                    result['formulas'][formula.name] = {
                        'value': bcfw_float,
                        'ratio': None
                    }
            except Exception:
                result['formulas'][formula.name] = None
        
    except Exception:
        pass
    
    return result


def run_phase2_parallel(n_trials=10000, n_cores=8):
    """
    Run Phase 2 tests in parallel.
    
    Returns dict with results for all formulas.
    """
    log_progress("PHASE2", f"Starting Phase 2: BCFW Amplituhedron ({n_trials} trials, {n_cores} cores)")
    
    # Initialize formulas
    formulas = [
        StandardBCFW(),
        GravityDoubledPT(),
        KLTBilinear()
    ]
    
    log_progress("PHASE2", f"Testing {len(formulas)} formula variants")
    
    # Check for checkpoint
    checkpoint = load_checkpoint('phase2_checkpoint.json')
    start_seed = 0
    results = []
    
    if checkpoint:
        log_progress("PHASE2", f"Resuming from checkpoint: {len(checkpoint.get('results', []))} points already tested")
        results = checkpoint.get('results', [])
        start_seed = checkpoint.get('last_seed', 0) + 1
    
    # Generate seed range
    seeds = list(range(start_seed, start_seed + n_trials - len(results)))
    
    log_progress("PHASE2", f"Testing {len(seeds)} new points...")
    
    # Run in parallel
    t_start = time.time()
    with Pool(processes=n_cores) as pool:
        new_results = pool.starmap(test_single_point_bcfw, [(seed, formulas) for seed in seeds])
    
    results.extend(new_results)
    elapsed = time.time() - t_start
    
    log_progress("PHASE2", f"Completed {len(results)} tests in {elapsed:.1f}s")
    
    # Save checkpoint
    save_checkpoint({
        'results': results,
        'last_seed': seeds[-1] if seeds else start_seed - 1,
        'timestamp': time.time()
    }, 'phase2_checkpoint.json')
    
    # Analyze results
    successful = [r for r in results if r.get('hodges') is not None]
    log_progress("PHASE2", f"Successful tests: {len(successful)}/{len(results)}")
    
    if len(successful) < 100:
        log_progress("PHASE2", "WARNING: Too few successful tests", "WARNING")
        return {
            'success': False,
            'error': 'insufficient_data',
            'successful_count': len(successful)
        }
    
    # Analyze each formula
    formula_results = {}
    
    for formula in formulas:
        formula_name = formula.name
        ratios = []
        exact_matches = 0
        
        for r in successful:
            if formula_name in r.get('formulas', {}):
                formula_data = r['formulas'][formula_name]
                if formula_data and formula_data.get('ratio') is not None:
                    ratio = formula_data['ratio']
                    ratios.append(ratio)
                    if abs(ratio - 1.0) < 1e-10:
                        exact_matches += 1
        
        if ratios:
            stats = statistical_analysis(ratios, formula_name)
            is_const, const_stats = is_constant_ratio(ratios, threshold=0.01)
            
            formula_results[formula_name] = {
                'stats': stats,
                'is_constant': is_const,
                'exact_matches': exact_matches,
                'count': len(ratios)
            }
            
            log_progress("PHASE2", f"{formula_name}: {len(ratios)} valid, {exact_matches} exact matches, cv={stats['cv']:.6e}")
            
            if is_const:
                norm = 1.0 / stats['mean']
                log_progress("PHASE2", f"  → Constant ratio! Normalization: {norm:.6e}")
        else:
            formula_results[formula_name] = {
                'error': 'no_valid_results'
            }
    
    # Find best formula
    best_formula = None
    best_score = -1
    
    for formula_name, data in formula_results.items():
        if 'error' in data:
            continue
        
        # Score: prefer exact matches, then constant ratios
        score = 0
        if data.get('exact_matches', 0) > 0:
            score = 1000 + data['exact_matches']
        elif data.get('is_constant', False):
            score = 100 + data.get('count', 0)
        else:
            score = data.get('count', 0)
        
        if score > best_score:
            best_score = score
            best_formula = formula_name
    
    phase2_result = {
        'success': True,
        'formula_results': formula_results,
        'best_formula': best_formula,
        'best_score': best_score,
        'successful_count': len(successful),
        'total_tests': len(results)
    }
    
    log_progress("PHASE2", f"Phase 2 complete. Best formula: {best_formula}")
    
    return phase2_result


def main():
    """Run Phase 2 standalone."""
    n_cores = min(8, cpu_count())
    result = run_phase2_parallel(n_trials=10000, n_cores=n_cores)
    
    # Save final result
    save_checkpoint(result, 'phase2_result.json')
    
    print("\n" + "="*70)
    print("PHASE 2 RESULTS")
    print("="*70)
    print(f"Best formula: {result.get('best_formula', 'none')}")
    print(f"Successful tests: {result.get('successful_count', 0)}")
    print("="*70)
    
    return result


if __name__ == '__main__':
    main()

