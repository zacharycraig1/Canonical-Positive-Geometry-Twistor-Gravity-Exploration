#!/usr/bin/env sage
# =============================================================================
# MHV SINGLE SOLUTION ANALYSIS
# =============================================================================
# Investigates the Du-Teng-Wu theory (arXiv:1603.08158) that only ONE of the
# (n-3)!=6 solutions contributes for MHV amplitudes.
#
# This addresses the ~2220x CHY vs Hodges discrepancy by:
# 1. Verifying that only 1 solution has non-zero Pfaffian
# 2. Checking normalization conventions
# 3. Comparing single-solution CHY to Hodges

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/psi_matrix.sage")
load("src/gravity_proof/positivity_region.sage")
load("src/gravity_proof/scattering_solver.sage")
load("src/gravity_proof/canonical_form_gravity.sage")
load("src/gravity_proof/amplitude_comparison.sage")

# Import Python modules
from src.kinematics.spinors import SpinorKinematics


def main(seed=42, verbose=True):
    """
    Run complete MHV single-solution analysis.
    
    Args:
        seed: Random seed for kinematics
        verbose: If True, print detailed output
        
    Returns:
        dict with all analysis results
    """
    if verbose:
        print("="*80)
        print("MHV SINGLE SOLUTION ANALYSIS")
        print("Investigating CHY Normalization via Du-Teng-Wu arXiv:1603.08158")
        print("="*80)
    
    # ==========================================================================
    # SETUP
    # ==========================================================================
    if verbose:
        print("\n[Setup] Creating kinematics and solving scattering equations...")
    
    kin = SpinorKinematics.random_rational(6, seed=seed)
    solver = ScatteringEquationSolver(kin)
    solutions = solver.solve_numerical()
    
    if not solutions:
        print("ERROR: No solutions found!")
        return {'status': 'failed', 'error': 'No solutions found'}
    
    if verbose:
        print(f"Found {len(solutions)} solutions (expected 6 for n=6)")
    
    # Create Psi matrix and comparison object
    var('z4 z5 z6')
    psi = PsiMatrixMHV(kin, z4, z5, z6)
    region = PositivityRegionR6(psi)
    canon = CanonicalFormR6(region)
    comp = AmplitudeComparison(canon, kin)
    
    results = {
        'seed': seed,
        'num_solutions': len(solutions),
        'tasks': {}
    }
    
    # ==========================================================================
    # TASK 1: Analyze Per-Solution Contributions
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("TASK 1: ANALYZE PER-SOLUTION CONTRIBUTIONS")
        print("="*80)
        print("Compute Pf'(Ψ) at each solution to identify which contribute")
        print("-"*80)
    
    try:
        per_solution_results = comp.analyze_per_solution_contributions(solver)
        results['tasks']['per_solution_analysis'] = {
            'status': 'success',
            'num_nonzero': sum(1 for c in per_solution_results['contributions'] if abs(c['contribution']) > 1e-10),
            'dominant_solution_index': per_solution_results['dominant_solution_index'],
            'dominant_contribution': per_solution_results['dominant_contribution']
        }
    except Exception as e:
        if verbose:
            print(f"ERROR in per-solution analysis: {e}")
            import traceback
            traceback.print_exc()
        results['tasks']['per_solution_analysis'] = {'status': 'failed', 'error': str(e)}
        per_solution_results = None
    
    # ==========================================================================
    # TASK 2: Verify MHV Single Solution Theory
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("TASK 2: VERIFY DU-TENG-WU MHV SINGLE SOLUTION THEORY")
        print("="*80)
        print("Theory: Exactly 1 of 6 solutions should have Pf'(Ψ) ≠ 0")
        print("-"*80)
    
    try:
        mhv_theory_results = comp.verify_mhv_single_solution_theory(solver)
        results['tasks']['mhv_theory_verification'] = {
            'status': 'success',
            'theory_confirmed': mhv_theory_results.get('theory_confirmed', False),
            'num_zero_pf': mhv_theory_results.get('num_zero_pf', 0),
            'num_nonzero_pf': mhv_theory_results.get('num_nonzero_pf', 0),
            'mhv_solution_matches_hodges': mhv_theory_results.get('mhv_solution_matches_hodges', False)
        }
    except Exception as e:
        if verbose:
            print(f"ERROR in MHV theory verification: {e}")
            import traceback
            traceback.print_exc()
        results['tasks']['mhv_theory_verification'] = {'status': 'failed', 'error': str(e)}
        mhv_theory_results = None
    
    # ==========================================================================
    # TASK 3: Compute CHY (MHV Solution Only)
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("TASK 3: COMPUTE CHY AMPLITUDE (MHV SOLUTION ONLY)")
        print("="*80)
        print("Use only the dominant solution and compare to Hodges")
        print("-"*80)
    
    try:
        single_solution_results = comp.compute_chy_mhv_single_solution(solver)
        results['tasks']['single_solution_chy'] = {
            'status': 'success',
            'dominant_solution_index': single_solution_results.get('dominant_solution_index', -1),
            'dominant_contribution': single_solution_results.get('dominant_contribution', 0),
            'hodges_amplitude': single_solution_results.get('hodges_amplitude', 0),
            'match': single_solution_results.get('match', False),
            'relative_difference': single_solution_results.get('relative_difference', float('inf'))
        }
    except Exception as e:
        if verbose:
            print(f"ERROR in single solution CHY: {e}")
            import traceback
            traceback.print_exc()
        results['tasks']['single_solution_chy'] = {'status': 'failed', 'error': str(e)}
        single_solution_results = None
    
    # ==========================================================================
    # TASK 4: Check Reduced Pfaffian Prefactor
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("TASK 4: CHECK REDUCED PFAFFIAN PREFACTOR")
        print("="*80)
        print("Verify sign factor (-1)^(i+j) and denominator 1/(z_i - z_j)")
        print("-"*80)
    
    try:
        # Pick the dominant solution for detailed check
        if per_solution_results and 'sorted_contributions' in per_solution_results:
            dominant_idx = per_solution_results['sorted_contributions'][0]['index']
            sol = solutions[dominant_idx]
            z4_val, z5_val, z6_val = sol['z4'], sol['z5'], sol['z6']
            
            psi_at_sol = PsiMatrixMHV(kin, z4_val, z5_val, z6_val)
            
            # Test different deletion conventions
            deletions_to_test = [
                (0, 1, "particles 1,2"),
                (0, 5, "particles 1,6"),
                (2, 5, "particles 3,6")
            ]
            
            pfaffian_results = {}
            if verbose:
                print(f"\nTesting deletion conventions on solution {dominant_idx+1}:")
            
            for del_i, del_j, description in deletions_to_test:
                try:
                    if del_i == 2 or del_j == 2:
                        # Use special method for deletion involving particle 3 (z3=∞)
                        pf = psi_at_sol.reduced_pfaffian_delete_3_6()
                    else:
                        pf = psi_at_sol.reduced_pfaffian_standard((del_i, del_j))
                    
                    try:
                        pf_float = float(pf) if hasattr(pf, '__float__') else complex(pf).real
                    except:
                        pf_float = pf
                    
                    pfaffian_results[f"delete_{del_i}_{del_j}"] = pf_float
                    if verbose:
                        print(f"  Delete {description}: Pf' = {pf_float:.6e}")
                except Exception as e:
                    pfaffian_results[f"delete_{del_i}_{del_j}"] = None
                    if verbose:
                        print(f"  Delete {description}: FAILED - {e}")
            
            results['tasks']['pfaffian_prefactor_check'] = {
                'status': 'success',
                'tested_deletions': pfaffian_results
            }
        else:
            results['tasks']['pfaffian_prefactor_check'] = {
                'status': 'skipped',
                'reason': 'No dominant solution identified'
            }
    except Exception as e:
        if verbose:
            print(f"ERROR in Pfaffian prefactor check: {e}")
            import traceback
            traceback.print_exc()
        results['tasks']['pfaffian_prefactor_check'] = {'status': 'failed', 'error': str(e)}
    
    # ==========================================================================
    # TASK 5: Check det'(Φ) Denominator Convention
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("TASK 5: CHECK det'(Φ) DENOMINATOR CONVENTION")
        print("="*80)
        print("Formula: det'(Φ) = det(Φ^{ijk}_{ijk}) / (σ_ij σ_jk σ_ki)²")
        print("-"*80)
    
    try:
        if solutions:
            sol = solutions[0]
            z4_val, z5_val, z6_val = sol['z4'], sol['z5'], sol['z6']
            
            # Compute det using solver method
            det_solver = solver.detprime_phi(z4_val, z5_val, z6_val)
            
            # Compute raw Jacobian determinant
            J = solver.jacobian_matrix(z4_val, z5_val, z6_val)
            det_raw = J.det()
            
            # Check if Vandermonde factor is included
            det_raw_alt, vander, det_prime_alt = solver.detprime_phi_with_vandermonde(z4_val, z5_val, z6_val)
            
            try:
                det_solver_float = float(det_solver) if hasattr(det_solver, '__float__') else complex(det_solver).real
                det_raw_float = float(det_raw) if hasattr(det_raw, '__float__') else complex(det_raw).real
                det_prime_alt_float = float(det_prime_alt) if hasattr(det_prime_alt, '__float__') else complex(det_prime_alt).real
            except:
                det_solver_float = det_solver
                det_raw_float = det_raw
                det_prime_alt_float = det_prime_alt
            
            if verbose:
                print(f"\nSolution 1:")
                print(f"  det'(Φ) via solver.detprime_phi(): {det_solver_float:.6e}")
                print(f"  det(Jacobian) raw:                 {det_raw_float:.6e}")
                print(f"  Ratio:                             {float(det_solver_float/det_raw_float) if det_raw_float != 0 else 'inf':.6f}")
                print(f"\nAlternative computation:")
                print(f"  det(Jacobian):        {det_raw_alt:.6e}")
                print(f"  Vandermonde factor:   {vander}")
                print(f"  det' (with Vander):   {det_prime_alt_float:.6e}")
            
            results['tasks']['detprime_convention_check'] = {
                'status': 'success',
                'det_solver': det_solver_float,
                'det_raw': det_raw_float,
                'det_prime_alt': det_prime_alt_float,
                'vandermonde_factor': vander,
                'ratio_solver_raw': float(det_solver_float/det_raw_float) if det_raw_float != 0 else None
            }
        else:
            results['tasks']['detprime_convention_check'] = {
                'status': 'skipped',
                'reason': 'No solutions available'
            }
    except Exception as e:
        if verbose:
            print(f"ERROR in det' convention check: {e}")
            import traceback
            traceback.print_exc()
        results['tasks']['detprime_convention_check'] = {'status': 'failed', 'error': str(e)}
    
    # ==========================================================================
    # TASK 6: Test at n=4 (Simpler Case)
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("TASK 6: TEST AT n=4 (ONLY 1 SOLUTION)")
        print("="*80)
        print("Simpler test case to isolate convention issues")
        print("-"*80)
        print("Note: This requires n=4 implementation - SKIPPED for now")
    
    results['tasks']['n4_test'] = {
        'status': 'skipped',
        'reason': 'n=4 implementation not available in current framework'
    }
    
    # ==========================================================================
    # SUMMARY AND CONCLUSIONS
    # ==========================================================================
    if verbose:
        print("\n" + "="*80)
        print("SUMMARY AND CONCLUSIONS")
        print("="*80)
    
    # Compile summary
    summary = {
        'mhv_theory_confirmed': False,
        'mhv_solution_matches_hodges': False,
        'normalization_factor_identified': None,
        'recommendations': []
    }
    
    # Check MHV theory
    if mhv_theory_results:
        summary['mhv_theory_confirmed'] = mhv_theory_results.get('theory_confirmed', False)
        summary['mhv_solution_matches_hodges'] = mhv_theory_results.get('mhv_solution_matches_hodges', False)
    
    # Print conclusions
    if verbose:
        if summary['mhv_theory_confirmed']:
            print("\n✅ MHV SINGLE SOLUTION THEORY CONFIRMED")
            print(f"   Only {mhv_theory_results.get('num_nonzero_pf', 0)} of {len(solutions)} solutions has non-zero Pf'(Ψ)")
            
            if summary['mhv_solution_matches_hodges']:
                print("   ✅ MHV solution contribution matches Hodges amplitude!")
                print("\n   CONCLUSION: CHY normalization issue is RESOLVED")
                summary['recommendations'].append("Use single MHV solution for amplitude calculation")
            else:
                print("   ⚠️  MHV solution does NOT match Hodges amplitude")
                if single_solution_results:
                    ratio = single_solution_results.get('ratio', 0)
                    print(f"   Ratio (single/Hodges) = {ratio:.6f}")
                print("\n   NEXT STEPS:")
                print("   - Verify reduced Pfaffian sign convention")
                print("   - Check det'(Φ) normalization factor")
                print("   - Test alternative gauge fixings")
                summary['recommendations'].extend([
                    "Check Pfaffian sign convention",
                    "Verify det'(Φ) normalization",
                    "Test alternative gauge fixings"
                ])
        else:
            print("\n⚠️  MHV SINGLE SOLUTION THEORY NOT CONFIRMED")
            if mhv_theory_results:
                print(f"   Found {mhv_theory_results.get('num_nonzero_pf', 0)} solutions with non-zero Pf'")
                print(f"   (Expected: exactly 1)")
            print("\n   POSSIBLE CAUSES:")
            print("   - Numerical threshold too tight")
            print("   - Non-MHV kinematics (check helicity configuration)")
            print("   - Implementation error in Pf' computation")
            summary['recommendations'].extend([
                "Check numerical thresholds",
                "Verify MHV helicity configuration",
                "Debug Pf' computation"
            ])
    
    results['summary'] = summary
    results['status'] = 'success'
    
    # Final status message
    if verbose:
        print("\n" + "="*80)
        print("ANALYSIS COMPLETE")
        print("="*80)
    
    return results


def run_multiple_kinematics(num_tests=5, seed_start=0, verbose=False):
    """
    Run MHV single-solution analysis with multiple different kinematics.
    
    This helps determine if results are generic or specific to particular kinematics.
    
    Args:
        num_tests: Number of different kinematic configurations to test
        seed_start: Starting seed value
        verbose: Print detailed output for each test
        
    Returns:
        List of result dicts, one per test
    """
    print("="*80)
    print(f"MHV SINGLE SOLUTION ANALYSIS - MULTIPLE KINEMATICS")
    print(f"Testing {num_tests} different kinematic configurations")
    print("="*80)
    
    all_results = []
    success_count = 0
    mhv_theory_count = 0
    hodges_match_count = 0
    
    for i in range(num_tests):
        seed = seed_start + i
        print(f"\n{'='*80}")
        print(f"Test {i+1}/{num_tests} (seed={seed})")
        print(f"{'='*80}")
        
        result = main(seed=seed, verbose=verbose)
        all_results.append(result)
        
        if result.get('status') == 'success':
            success_count += 1
            
            # Check MHV theory confirmation
            if result.get('summary', {}).get('mhv_theory_confirmed', False):
                mhv_theory_count += 1
            
            # Check Hodges match
            if result.get('summary', {}).get('mhv_solution_matches_hodges', False):
                hodges_match_count += 1
    
    # Overall summary
    print("\n" + "="*80)
    print("OVERALL SUMMARY (ALL KINEMATICS)")
    print("="*80)
    print(f"Successful runs: {success_count}/{num_tests}")
    print(f"MHV theory confirmed: {mhv_theory_count}/{success_count}")
    print(f"MHV solution = Hodges: {hodges_match_count}/{success_count}")
    
    if hodges_match_count == success_count:
        print("\n✅ ALL TESTS: MHV solution matches Hodges!")
        print("   Normalization issue is RESOLVED")
    elif mhv_theory_count == success_count:
        print("\n⚠️  MHV theory confirmed in all tests, but normalization factor needed")
        print(f"   {hodges_match_count}/{success_count} tests match Hodges")
    else:
        print("\n❌ MHV theory NOT consistently confirmed")
        print(f"   Only {mhv_theory_count}/{success_count} tests confirm MHV theory")
    
    return all_results


if __name__ == "__main__":
    # Run analysis with single kinematics
    results = main(seed=42, verbose=True)
    
    # Uncomment to test multiple kinematics
    # all_results = run_multiple_kinematics(num_tests=3, seed_start=42, verbose=False)

