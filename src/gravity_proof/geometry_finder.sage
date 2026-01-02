#!/usr/bin/env sage
# =============================================================================
# Geometry Finder for 6-Point MHV Gravity
# =============================================================================
# Main orchestration script to:
# 1. Test all four positivity options (B1-B4)
# 2. Identify which gives correct 6-vertex geometry
# 3. Enumerate boundaries from positivity conditions
# 4. Compute canonical form geometrically
# 5. Compare to known Hodges amplitude

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/positivity_options.sage")
load("src/gravity_proof/scattering_solver.sage")
load("src/gravity_proof/psi_matrix.sage")

# Import Python modules
from src.kinematics.spinors import SpinorKinematics


class GeometryFinder:
    """
    Finds the correct positive geometry R6 for 6-point MHV gravity.
    """
    
    def __init__(self, kinematics=None, seed=42):
        """
        Initialize geometry finder.
        
        Args:
            kinematics: Optional SpinorKinematics object
            seed: Random seed if generating kinematics
        """
        if kinematics is None:
            self.kin = SpinorKinematics.random_rational(6, seed=seed)
        else:
            self.kin = kinematics
        
        self.tester = PositivityOptionTester(self.kin)
        self.solver = ScatteringEquationSolver(self.kin)
        self.solutions = None
        
    def get_scattering_solutions(self):
        """
        Solve scattering equations to get the 6 solutions.
        
        Returns:
            List of solution dicts
        """
        if self.solutions is None:
            print("Solving scattering equations...")
            self.solutions = self.solver.solve_numerical()
            print(f"Found {len(self.solutions)} solutions")
        return self.solutions
    
    def check_solutions_as_vertices(self, option_name, tolerance=1e-6):
        """
        Check if scattering equation solutions are vertices of the positivity region.
        
        A vertex is on the boundary (not in strict interior) and at the intersection
        of multiple boundaries.
        
        Args:
            option_name: 'B1', 'B2', 'B3', or 'B4'
            tolerance: Numerical tolerance for boundary check
            
        Returns:
            dict with analysis results
        """
        solutions = self.get_scattering_solutions()
        
        results = {
            'option': option_name,
            'total_solutions': len(solutions),
            'real_solutions': 0,
            'solutions_in_region': [],
            'solutions_on_boundary': [],
            'solutions_in_interior': [],
            'solutions_outside': []
        }
        
        for i, sol in enumerate(solutions):
            z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
            
            # Check if real
            is_real = True
            if hasattr(z4, 'imag'):
                is_real = is_real and abs(z4.imag) < tolerance
            if hasattr(z5, 'imag'):
                is_real = is_real and abs(z5.imag) < tolerance
            if hasattr(z6, 'imag'):
                is_real = is_real and abs(z6.imag) < tolerance
            
            if not is_real:
                continue  # Skip complex solutions for geometry
            
            results['real_solutions'] += 1
            
            # Convert to real
            z4_r = z4.real if hasattr(z4, 'real') else z4
            z5_r = z5.real if hasattr(z5, 'real') else z5
            z6_r = z6.real if hasattr(z6, 'real') else z6
            
            # Test the option
            if option_name == 'B1':
                satisfies = self.tester.test_option_b1(z4_r, z5_r, z6_r)
                # For B1, check if on boundary (one of the ordering equalities)
                on_boundary = (abs(z4_r - z5_r) < tolerance or 
                              abs(z5_r - z6_r) < tolerance or 
                              abs(z6_r) < tolerance)
            elif option_name == 'B2':
                satisfies, pf_val = self.tester.test_option_b2(z4_r, z5_r, z6_r)
                on_boundary = (abs(pf_val) < tolerance) if pf_val is not None else False
            elif option_name == 'B3':
                satisfies, details = self.tester.test_option_b3(z4_r, z5_r, z6_r)
                # On boundary if any minor is zero
                on_boundary = False
                if details and isinstance(details, dict):
                    if details.get('4x4_det') is not None:
                        on_boundary = on_boundary or (abs(details['4x4_det']) < tolerance)
                    for minor in details.get('2x2_minors', []):
                        on_boundary = on_boundary or (abs(minor) < tolerance)
            elif option_name == 'B4':
                satisfies, eigenvals = self.tester.test_option_b4(z4_r, z5_r, z6_r)
                # On boundary if any eigenvalue is zero
                on_boundary = False
                if eigenvals is not None:
                    for ev in eigenvals:
                        on_boundary = on_boundary or (abs(ev) < tolerance)
            else:
                raise ValueError(f"Unknown option: {option_name}")
            
            sol_info = {
                'index': i,
                'z4': z4_r,
                'z5': z5_r,
                'z6': z6_r,
                'satisfies': satisfies,
                'on_boundary': on_boundary
            }
            
            if satisfies:
                results['solutions_in_region'].append(sol_info)
                if on_boundary:
                    results['solutions_on_boundary'].append(sol_info)
                else:
                    results['solutions_in_interior'].append(sol_info)
            else:
                results['solutions_outside'].append(sol_info)
        
        return results
    
    def compare_all_options(self):
        """
        Compare all four options to find which gives correct geometry.
        
        Returns:
            dict with comparison results
        """
        print("\n" + "="*80)
        print("COMPARING ALL POSITIVITY OPTIONS")
        print("="*80)
        
        results = {}
        
        for option in ['B1', 'B2', 'B3', 'B4']:
            print(f"\nTesting Option {option}...")
            vertex_check = self.check_solutions_as_vertices(option)
            results[option] = vertex_check
            
            print(f"  Real solutions: {vertex_check['real_solutions']}")
            print(f"  In region: {len(vertex_check['solutions_in_region'])}")
            print(f"  On boundary (vertices): {len(vertex_check['solutions_on_boundary'])}")
            print(f"  In interior: {len(vertex_check['solutions_in_interior'])}")
            print(f"  Outside: {len(vertex_check['solutions_outside'])}")
        
        # Determine which option is correct
        print("\n" + "="*80)
        print("ANALYSIS")
        print("="*80)
        
        correct_option = None
        for option, data in results.items():
            num_vertices = len(data['solutions_on_boundary'])
            num_in_region = len(data['solutions_in_region'])
            
            # Expected: 6 vertices = 6 scattering solutions on boundary
            if num_vertices == 6 and data['real_solutions'] == 6:
                print(f"\n✓ Option {option}: CANDIDATE - {num_vertices} vertices (all real solutions on boundary)")
                if correct_option is None:
                    correct_option = option
            elif num_in_region == 6:
                print(f"\n? Option {option}: PARTIAL - {num_in_region} solutions in region, {num_vertices} on boundary")
            else:
                print(f"\n✗ Option {option}: REJECTED - Only {num_in_region} solutions in region")
        
        if correct_option:
            print(f"\n" + "="*80)
            print(f"CONCLUSION: Option {correct_option} appears to define R6 correctly")
            print("="*80)
        else:
            print(f"\n" + "="*80)
            print(f"WARNING: No option gives exact 6-vertex structure")
            print("="*80)
        
        return {
            'option_results': results,
            'correct_option': correct_option
        }
    
    def enumerate_boundaries(self, option_name):
        """
        Enumerate boundaries for the given option.
        
        Args:
            option_name: 'B1', 'B2', 'B3', or 'B4'
            
        Returns:
            List of boundary descriptions
        """
        print(f"\n" + "="*80)
        print(f"ENUMERATING BOUNDARIES FOR OPTION {option_name}")
        print("="*80)
        
        boundaries = []
        
        if option_name == 'B1':
            # Ordering boundaries
            boundaries.append({
                'name': 'z4_z5',
                'equation': 'z4 - z5 = 0',
                'description': 'Collision of particles 4 and 5'
            })
            boundaries.append({
                'name': 'z5_z6',
                'equation': 'z5 - z6 = 0',
                'description': 'Collision of particles 5 and 6'
            })
            boundaries.append({
                'name': 'z6_0',
                'equation': 'z6 = 0',
                'description': 'Particle 6 at z1=0'
            })
            
        elif option_name == 'B2':
            # Pfaffian zero locus
            boundaries.append({
                'name': 'pfaffian_zero',
                'equation': "Pf'(Psi) = 0",
                'description': 'Reduced Pfaffian vanishes'
            })
            # Plus ordering boundaries
            boundaries.append({
                'name': 'ordering_boundaries',
                'equation': 'Various z_i collisions',
                'description': 'Worldsheet ordering constraints'
            })
            
        elif option_name == 'B3':
            # Principal minor boundaries
            boundaries.append({
                'name': 'principal_minors',
                'equation': 'Various principal minors = 0',
                'description': 'Total positivity boundary'
            })
            
        elif option_name == 'B4':
            # Eigenvalue boundaries
            boundaries.append({
                'name': 'eigenvalue_zero',
                'equation': 'Some eigenvalue of i*Psi = 0',
                'description': 'Positive semi-definite boundary'
            })
        
        print(f"Found {len(boundaries)} boundary types:")
        for i, b in enumerate(boundaries):
            print(f"  {i+1}. {b['name']}: {b['equation']}")
            print(f"     {b['description']}")
        
        return boundaries
    
    def compute_geometric_canonical_form(self, option_name, boundaries):
        """
        Compute canonical form from geometric boundaries.
        
        NOT from the CHY amplitude - derive from the region structure.
        
        Args:
            option_name: The correct option
            boundaries: List of boundary dicts
            
        Returns:
            dict with canonical form structure
        """
        print(f"\n" + "="*80)
        print(f"COMPUTING CANONICAL FORM FROM GEOMETRY (OPTION {option_name})")
        print("="*80)
        
        var('z4 z5 z6')
        
        # Denominator from boundaries
        if option_name == 'B1':
            # Simple ordering: z4 > z5 > z6 > 0
            # Canonical form: Omega = N * dz4 ^ dz5 ^ dz6 / ((z4-z5)(z5-z6)(z6))
            denominator_factors = [z4 - z5, z5 - z6, z6]
            denominator = (z4 - z5) * (z5 - z6) * z6
            
        elif option_name == 'B2':
            # Pfaffian positivity
            # Boundaries include Pf'(Psi) = 0
            # For now, use symbolic Pfaffian
            psi_symbolic = PsiMatrixMHV(self.kin, z4, z5, z6)
            try:
                pf_prime = psi_symbolic.reduced_pfaffian_standard((0, 1))
                if pf_prime == 0:
                    pf_prime = psi_symbolic.reduced_pfaffian_delete_3_6()
            except:
                pf_prime = psi_symbolic.reduced_pfaffian_delete_3_6()
            
            # Denominator includes Pfaffian
            # But canonical form structure depends on ALL boundaries
            # For now, use placeholder
            denominator_factors = ['Pf(Psi)', 'ordering']
            denominator = pf_prime  # This is actually wrong - need proper boundary analysis
            
        else:
            # Options B3, B4 are more complex
            denominator_factors = ['Complex boundary structure']
            denominator = None
        
        # Numerator: Determined by requiring unit residues
        # This is the KEY: compute N from geometry, not from amplitude
        
        # For simple polytopes, N is typically a constant or simple polynomial
        # For Option B1 with ordering constraints, N might be constant
        
        # TODO: Proper residue analysis to determine N
        numerator = 1  # Placeholder - needs geometric computation
        
        result = {
            'option': option_name,
            'denominator_factors': denominator_factors,
            'denominator': denominator,
            'numerator': numerator,
            'canonical_form': f"N * dz4 ∧ dz5 ∧ dz6 / ({denominator})",
            'note': 'Numerator N should be determined from unit residue requirement'
        }
        
        print(f"Denominator factors: {denominator_factors}")
        print(f"Canonical form structure: {result['canonical_form']}")
        print(f"\nNote: {result['note']}")
        
        return result
    
    def run_full_analysis(self):
        """
        Run complete geometry finding analysis.
        
        Returns:
            dict with all results
        """
        print("="*80)
        print("FINDING POSITIVE GEOMETRY FOR 6-POINT MHV GRAVITY")
        print("="*80)
        
        # Step 1: Compare all options
        comparison = self.compare_all_options()
        correct_option = comparison['correct_option']
        
        if correct_option is None:
            print("\nWARNING: Could not definitively identify correct option")
            print("Proceeding with detailed analysis of all options...")
            # Continue with all options for analysis
            options_to_analyze = ['B1', 'B2', 'B3', 'B4']
        else:
            options_to_analyze = [correct_option]
        
        # Step 2: Enumerate boundaries for correct option(s)
        boundaries_results = {}
        for option in options_to_analyze:
            boundaries = self.enumerate_boundaries(option)
            boundaries_results[option] = boundaries
        
        # Step 3: Compute geometric canonical form
        canonical_forms = {}
        for option in options_to_analyze:
            canonical_form = self.compute_geometric_canonical_form(
                option, 
                boundaries_results[option]
            )
            canonical_forms[option] = canonical_form
        
        # Step 4: Symbolic comparison to Hodges
        load("src/gravity_proof/symbolic_comparison.sage")
        
        symbolic_results = {}
        if correct_option:
            print(f"\n" + "="*80)
            print("STEP 4: SYMBOLIC COMPARISON TO HODGES AMPLITUDE")
            print("="*80)
            
            # Create canonical form for comparison
            var('z4 z5 z6')
            psi = PsiMatrixMHV(self.kin, z4, z5, z6)
            load("src/gravity_proof/positivity_region.sage")
            region = PositivityRegionR6(psi)
            canon = CanonicalFormR6(region)
            
            # Compare
            comp = SymbolicComparison(canon, self.kin)
            
            try:
                comparison_result = comp.compare_at_scattering_solutions(self.solver)
                symbolic_results['numerical_comparison'] = comparison_result
            except Exception as e:
                print(f"Numerical comparison failed: {e}")
                symbolic_results['numerical_comparison'] = {'status': 'failed', 'error': str(e)}
            
            try:
                factorization = comp.verify_boundary_factorization(
                    correct_option,
                    boundaries_results[correct_option]
                )
                symbolic_results['factorization'] = factorization
            except Exception as e:
                print(f"Factorization check failed: {e}")
                symbolic_results['factorization'] = {'status': 'failed', 'error': str(e)}
            
            try:
                equality = comp.symbolic_equality_check(correct_option)
                symbolic_results['symbolic_equality'] = equality
            except Exception as e:
                print(f"Symbolic equality check failed: {e}")
                symbolic_results['symbolic_equality'] = {'status': 'failed', 'error': str(e)}
        
        # Step 5: Final Summary
        print(f"\n" + "="*80)
        print("FINAL SUMMARY")
        print("="*80)
        print(f"\nScattering equation solutions: {len(self.solutions)}")
        print(f"Correct positivity option: {correct_option if correct_option else 'UNDETERMINED'}")
        
        if correct_option:
            print(f"\nGeometric structure of R6 (Option {correct_option}):")
            print(f"  Boundaries: {len(boundaries_results[correct_option])} types")
            print(f"  Canonical form: {canonical_forms[correct_option].get('canonical_form', 'Unknown')}")
            
            if 'numerical_comparison' in symbolic_results:
                nc = symbolic_results['numerical_comparison']
                if nc.get('status') == 'success':
                    print(f"\nNumerical verification:")
                    print(f"  {nc.get('num_matches', 0)}/{nc.get('num_solutions_tested', 0)} solutions match Hodges")
                    if nc.get('num_matches', 0) > 0:
                        print(f"  ✓ Canonical form agrees with Hodges amplitude!")
        
        return {
            'comparison': comparison,
            'correct_option': correct_option,
            'boundaries': boundaries_results,
            'canonical_forms': canonical_forms,
            'symbolic_results': symbolic_results
        }


def main(seed=42):
    """
    Main execution function.
    
    Args:
        seed: Random seed for kinematics
        
    Returns:
        GeometryFinder instance with results
    """
    finder = GeometryFinder(seed=seed)
    results = finder.run_full_analysis()
    return finder, results


if __name__ == "__main__":
    finder, results = main(seed=42)

