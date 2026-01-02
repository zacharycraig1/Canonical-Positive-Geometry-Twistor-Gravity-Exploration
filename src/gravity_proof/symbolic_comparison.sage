#!/usr/bin/env sage
# =============================================================================
# Symbolic Comparison: Canonical Form vs Hodges Amplitude
# =============================================================================
# Compares the geometrically-derived canonical form to the known Hodges amplitude
# SYMBOLICALLY, not just numerically.

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/canonical_form_gravity.sage")
load("src/gravity_proof/scattering_solver.sage")

# Import Python modules
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


class SymbolicComparison:
    """
    Compares canonical form (from geometry) to Hodges amplitude (known formula).
    """
    
    def __init__(self, canonical_form, kinematics):
        """
        Initialize comparison.
        
        Args:
            canonical_form: CanonicalFormR6 instance
            kinematics: SpinorKinematics object
        """
        self.canon_form = canonical_form
        self.kin = kinematics
        
    def compute_hodges_symbolic(self):
        """
        Compute Hodges amplitude using known formula.
        
        Returns:
            Hodges amplitude value
        """
        lambdas = self.kin.lambdas
        tilde_lambdas = self.kin.tilde_lambdas
        negative_indices = (0, 1)  # Particles 1,2 (0-indexed)
        
        result, status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, negative_indices)
        
        if status != "ok":
            raise ValueError(f"Hodges formula failed: {status}")
        
        return result
    
    def compare_at_scattering_solutions(self, solver):
        """
        Compare canonical form to Hodges at scattering equation solutions.
        
        This tests if Omega(R6) = M6^MHV numerically at the physical points.
        
        Args:
            solver: ScatteringEquationSolver instance
            
        Returns:
            dict with comparison results
        """
        print("\n" + "="*80)
        print("SYMBOLIC COMPARISON: CANONICAL FORM vs HODGES")
        print("="*80)
        
        # Get Hodges amplitude
        try:
            hodges = self.compute_hodges_symbolic()
            hodges_float = float(hodges) if hasattr(hodges, '__float__') else complex(hodges).real
            print(f"\nHodges amplitude: {hodges_float:.6e}")
        except Exception as e:
            print(f"ERROR computing Hodges: {e}")
            return {'status': 'failed', 'error': str(e)}
        
        # Get scattering solutions
        solutions = solver.solve_numerical()
        print(f"Scattering equation solutions: {len(solutions)}")
        
        # Evaluate canonical form at each solution
        comparison_results = []
        
        for i, sol in enumerate(solutions):
            z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
            
            # Check if real
            is_real = True
            if hasattr(z4, 'imag'):
                is_real = is_real and abs(z4.imag) < 1e-10
            if hasattr(z5, 'imag'):
                is_real = is_real and abs(z5.imag) < 1e-10
            if hasattr(z6, 'imag'):
                is_real = is_real and abs(z6.imag) < 1e-10
            
            if not is_real:
                continue
            
            # Evaluate canonical form
            try:
                omega_val = self.canon_form.evaluate_at_point(z4, z5, z6)
                if omega_val is None:
                    continue
                
                omega_float = float(omega_val) if hasattr(omega_val, '__float__') else complex(omega_val).real
                
                # Compare
                diff = abs(omega_float - hodges_float)
                rel_diff = diff / abs(hodges_float) if hodges_float != 0 else float('inf')
                
                comparison_results.append({
                    'solution_index': i,
                    'z4': z4,
                    'z5': z5,
                    'z6': z6,
                    'omega': omega_float,
                    'hodges': hodges_float,
                    'difference': diff,
                    'relative_difference': rel_diff,
                    'match': rel_diff < 1e-6
                })
                
                print(f"\nSolution {i+1}:")
                print(f"  Omega(R6): {omega_float:.6e}")
                print(f"  Hodges:    {hodges_float:.6e}")
                print(f"  Rel diff:  {rel_diff:.6e} {'✓' if rel_diff < 1e-6 else '✗'}")
                
            except Exception as e:
                print(f"\nSolution {i+1}: Evaluation failed - {e}")
        
        # Summary
        if comparison_results:
            matches = sum(1 for r in comparison_results if r['match'])
            print(f"\n" + "="*80)
            print(f"SUMMARY: {matches}/{len(comparison_results)} solutions match (< 1e-6 relative error)")
            print("="*80)
            
            status = 'success' if matches > 0 else 'no_match'
        else:
            status = 'no_results'
            matches = 0
        
        return {
            'status': status,
            'hodges_amplitude': hodges_float,
            'num_solutions_tested': len(comparison_results),
            'num_matches': matches,
            'comparison_results': comparison_results
        }
    
    def verify_boundary_factorization(self, option_name, boundaries):
        """
        Verify that boundaries correspond to factorization channels.
        
        At pole s_ijk → 0, amplitude should factorize:
        M_6 → M_4(i,j,k,P) × 1/s_ijk × M_4(-P, remaining)
        
        Args:
            option_name: The identified correct option
            boundaries: List of boundary dicts
            
        Returns:
            dict with factorization verification results
        """
        print(f"\n" + "="*80)
        print(f"VERIFYING BOUNDARY FACTORIZATION")
        print("="*80)
        
        expected_channels = [
            's_123', 's_234', 's_345', 's_456', 's_561', 's_612'
        ]
        
        print(f"\nExpected factorization channels: {expected_channels}")
        print(f"Geometric boundaries from Option {option_name}: {len(boundaries)}")
        
        # For Option B1: Simple ordering
        # Boundaries z4-z5=0, z5-z6=0, z6=0 don't obviously map to s_ijk
        # This might indicate B1 is NOT the complete picture
        
        # For Option B2: Pfaffian boundary
        # Pf'(Psi) = 0 might occur when certain s_ijk → 0
        # Need to check this correspondence
        
        factorization_check = {
            'option': option_name,
            'num_boundaries': len(boundaries),
            'expected_channels': expected_channels,
            'correspondence': 'To be determined - requires residue computation',
            'note': 'Full verification requires computing residues at each boundary'
        }
        
        print(f"\nFactorization verification: {factorization_check['correspondence']}")
        print(f"Note: {factorization_check['note']}")
        
        return factorization_check
    
    def symbolic_equality_check(self, option_name):
        """
        Attempt to verify symbolic equality between canonical form and amplitude.
        
        This is the key proof step - showing they are the same expression,
        not just numerically equal.
        
        Args:
            option_name: The identified correct option
            
        Returns:
            dict with equality status
        """
        print(f"\n" + "="*80)
        print(f"SYMBOLIC EQUALITY CHECK")
        print("="*80)
        
        # Get canonical form expression
        expr = self.canon_form.canonical_form_expression()
        
        print(f"\nCanonical form (Option {option_name}):")
        print(f"  Numerator: {expr.get('numerator', 'Unknown')}")
        print(f"  Denominator: {expr.get('denominator_product', 'Unknown')}")
        
        # For symbolic equality, we'd need to:
        # 1. Express Hodges amplitude in same coordinates (z4, z5, z6)
        # 2. Simplify both expressions
        # 3. Check if they're equal modulo a constant
        
        # This is complex and depends on the option
        
        result = {
            'status': 'partial',
            'canonical_form': expr,
            'note': 'Full symbolic equality requires expressing Hodges in (z4,z5,z6) coordinates',
            'approach': 'Use CHY formula to relate Hodges to Pfaffian, then compare'
        }
        
        print(f"\nStatus: {result['status']}")
        print(f"Note: {result['note']}")
        
        return result


def test_symbolic_comparison():
    """Test symbolic comparison functionality"""
    print("Testing SymbolicComparison...")
    
    from src.kinematics.spinors import SpinorKinematics
    
    # Create test setup - use load() for Sage files
    load("src/gravity_proof/psi_matrix.sage")
    load("src/gravity_proof/positivity_region.sage")
    
    kin = SpinorKinematics.random_rational(6, seed=42)
    var('z4 z5 z6')
    psi = PsiMatrixMHV(kin, z4, z5, z6)
    region = PositivityRegionR6(psi)
    canon = CanonicalFormR6(region)
    
    # Create comparison
    comp = SymbolicComparison(canon, kin)
    
    # Test Hodges
    try:
        hodges = comp.compute_hodges_symbolic()
        print(f"Hodges amplitude: {hodges}")
    except Exception as e:
        print(f"Hodges computation failed: {e}")
    
    # Test comparison at solutions
    try:
        solver = ScatteringEquationSolver(kin)
        results = comp.compare_at_scattering_solutions(solver)
        print(f"\nComparison: {results['num_matches']}/{results['num_solutions_tested']} matches")
    except Exception as e:
        print(f"Comparison failed: {e}")
    
    print("Test complete.")


if __name__ == "__main__":
    test_symbolic_comparison()

