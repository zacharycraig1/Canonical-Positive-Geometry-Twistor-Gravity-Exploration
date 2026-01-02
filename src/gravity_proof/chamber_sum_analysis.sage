#!/usr/bin/env sage
# =============================================================================
# Chamber Summation Analysis for 6-Point MHV Gravity
# =============================================================================
# Tests the hypothesis that gravity amplitude = SUM of canonical forms over
# all ordering chambers, not a single canonical form.
#
# Key insight: CHY sums over ALL (n-3)! = 6 solutions. Each solution belongs
# to a distinct ordering chamber. The "positive geometry" for gravity may be
# the UNION of all chambers.

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/scattering_solver.sage")
load("src/gravity_proof/psi_matrix.sage")

# Import Python modules
from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


class ChamberSumAnalysis:
    """
    Analyzes the chamber structure of scattering equation solutions
    and tests if the gravity amplitude is a sum over all chambers.
    
    The CHY formula for MHV gravity is:
        M_n = <12>^8 * Σ [Pf'(Ψ)]² / det'(Φ)
    
    where the helicity factor <12>^8 comes from the MHV configuration
    with particles 1,2 having negative helicity.
    """
    
    def __init__(self, kinematics=None, seed=42):
        """
        Initialize chamber analysis.
        
        Args:
            kinematics: Optional SpinorKinematics object
            seed: Random seed for kinematics generation
        """
        if kinematics is None:
            self.kin = SpinorKinematics.random_rational(6, seed=seed)
        else:
            self.kin = kinematics
        
        self.solver = ScatteringEquationSolver(self.kin)
        self.solutions = None
        self.chambers = None
        
        # Helicity factor <12>^8 for MHV with particles 1,2 negative helicity
        self.ang_12 = self.kin.angle(0, 1)
        self.helicity_factor = self.ang_12**8
    
    def solve_and_group_by_chamber(self):
        """
        Solve scattering equations and group solutions by ordering chamber.
        
        Each chamber is defined by an ordering of all 6 points on CP^1.
        With gauge fixing (z1, z2, z3) = (0, 1, infinity), we track the
        relative positions of z4, z5, z6 with respect to 0 and 1.
        
        Returns:
            dict mapping chamber ordering tuples to solution lists
        """
        print("="*80)
        print("CHAMBER SUMMATION ANALYSIS")
        print("="*80)
        
        # Solve scattering equations
        print("\nSolving scattering equations...")
        self.solutions = self.solver.solve_numerical()
        print(f"Found {len(self.solutions)} solutions")
        
        # Group by chamber
        self.chambers = {}
        
        for i, sol in enumerate(self.solutions):
            z4 = sol['z4']
            z5 = sol['z5']
            z6 = sol['z6']
            
            # Get real parts
            z4_r = self._to_real(z4)
            z5_r = self._to_real(z5)
            z6_r = self._to_real(z6)
            
            # Create ordering tuple including fixed points
            # z1 = 0, z2 = 1, z3 = infinity
            points = [
                (0, 'z1'),
                (1, 'z2'),
                (z4_r, 'z4'),
                (z5_r, 'z5'),
                (z6_r, 'z6')
            ]
            points_sorted = sorted(points, key=lambda x: x[0])
            ordering = tuple(p[1] for p in points_sorted)
            
            # Also store simple ordering of z4, z5, z6 only
            z_points = [(z4_r, 'z4'), (z5_r, 'z5'), (z6_r, 'z6')]
            z_sorted = sorted(z_points, key=lambda x: x[0], reverse=True)
            simple_ordering = ' > '.join([p[1] for p in z_sorted])
            
            if ordering not in self.chambers:
                self.chambers[ordering] = []
            
            self.chambers[ordering].append({
                'index': i,
                'z4': z4_r,
                'z5': z5_r,
                'z6': z6_r,
                'solution': sol,
                'simple_ordering': simple_ordering
            })
        
        # Print chamber structure
        print(f"\nFound {len(self.chambers)} distinct chambers:")
        for ordering, sols in sorted(self.chambers.items(), key=lambda x: str(x[0])):
            ordering_str = ' < '.join(ordering) + ' < z3(inf)'
            print(f"  {ordering_str}: {len(sols)} solution(s)")
            for s in sols:
                print(f"    Solution {s['index']+1}: z4={s['z4']:.4f}, z5={s['z5']:.4f}, z6={s['z6']:.4f}")
                print(f"      Simple ordering: {s['simple_ordering']}")
        
        return self.chambers
    
    def _to_real(self, z):
        """Convert to real float."""
        if hasattr(z, 'real') and callable(z.real):
            return float(z.real())
        elif hasattr(z, 'real'):
            return float(z.real)
        else:
            return float(z)
    
    def compute_per_chamber_contribution(self):
        """
        Compute the CHY contribution from each chamber.
        
        For each chamber, the contribution is:
            Omega_chamber = Pf'(Psi)^2 / det'(Phi)
        
        Returns:
            dict with per-chamber contributions
        """
        if self.chambers is None:
            self.solve_and_group_by_chamber()
        
        print("\n" + "="*80)
        print("PER-CHAMBER CONTRIBUTIONS")
        print("="*80)
        
        contributions = {}
        total_sum = 0
        
        for ordering, sols in self.chambers.items():
            ordering_str = ' < '.join(ordering)
            
            # Use first solution in chamber (should be only one per chamber)
            sol_data = sols[0]
            sol = sol_data['solution']
            z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
            
            # Compute Pf'(Psi)
            psi = PsiMatrixMHV(self.kin, z4, z5, z6)
            try:
                pf_prime = psi.reduced_pfaffian_standard((0, 1))
                if pf_prime == 0:
                    pf_prime = psi.reduced_pfaffian_delete_3_6()
            except:
                try:
                    pf_prime = psi.reduced_pfaffian_delete_3_6()
                except Exception as e:
                    print(f"  Chamber {ordering_str}: Pfaffian computation failed - {e}")
                    continue
            
            # Compute det'(Phi) - Jacobian determinant
            J = self.solver.jacobian_matrix(z4, z5, z6)
            det_J = J.det()
            
            if det_J == 0:
                print(f"  Chamber {ordering_str}: Zero Jacobian determinant")
                continue
            
            # CHY contribution: Pf'^2 / det(J)
            contribution = (pf_prime**2) / det_J
            
            # Convert to float for display
            try:
                contrib_float = float(contribution)
            except:
                contrib_float = complex(contribution).real
            
            contributions[ordering] = {
                'contribution': contribution,
                'contribution_float': contrib_float,
                'pf_prime': pf_prime,
                'det_J': det_J,
                'z4': z4,
                'z5': z5,
                'z6': z6
            }
            
            total_sum += contribution
            
            print(f"\n  Chamber: {ordering_str}")
            print(f"    Pf'(Psi) = {float(pf_prime):.6e}")
            print(f"    det'(Phi) = {float(det_J):.6e}")
            print(f"    Contribution = {contrib_float:.6e}")
        
        # Apply helicity factor <12>^8
        total_sum_with_helicity = self.helicity_factor * total_sum
        
        # Convert total to float
        try:
            total_float = float(total_sum_with_helicity)
        except:
            total_float = complex(total_sum_with_helicity).real
        
        print(f"\n" + "-"*40)
        print(f"SUM OVER ALL CHAMBERS (before helicity): {float(total_sum):.6e}")
        print(f"HELICITY FACTOR <12>^8: {float(self.helicity_factor):.6e}")
        print(f"TOTAL SUM (with helicity): {total_float:.6e}")
        
        return {
            'per_chamber': contributions,
            'total_sum': total_sum_with_helicity,
            'total_float': total_float,
            'num_chambers': len(contributions)
        }
    
    def compute_hodges_amplitude(self):
        """
        Compute the Hodges amplitude for comparison.
        
        Returns:
            dict with Hodges amplitude
        """
        print("\n" + "="*80)
        print("HODGES AMPLITUDE")
        print("="*80)
        
        lambdas = self.kin.lambdas
        tilde_lambdas = self.kin.tilde_lambdas
        negative_indices = (0, 1)  # Particles 1,2 have negative helicity
        
        try:
            result, status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, negative_indices)
            
            if status != "ok":
                print(f"Hodges computation failed: {status}")
                return {'status': 'failed', 'error': status}
            
            try:
                hodges_float = float(result)
            except:
                hodges_float = complex(result).real
            
            print(f"Hodges amplitude: {hodges_float:.6e}")
            
            return {
                'amplitude': result,
                'amplitude_float': hodges_float,
                'status': 'ok'
            }
        except Exception as e:
            print(f"Hodges computation error: {e}")
            return {'status': 'failed', 'error': str(e)}
    
    def verify_chamber_sum_equals_hodges(self):
        """
        Main verification: Does sum over chambers equal Hodges amplitude?
        
        Returns:
            dict with verification results
        """
        print("\n" + "="*80)
        print("VERIFICATION: CHAMBER SUM vs HODGES")
        print("="*80)
        
        # Compute chamber contributions
        chamber_result = self.compute_per_chamber_contribution()
        
        # Compute Hodges
        hodges_result = self.compute_hodges_amplitude()
        
        if hodges_result['status'] != 'ok':
            return {
                'status': 'failed',
                'error': 'Hodges computation failed',
                'chamber_sum': chamber_result['total_float'],
                'hodges': None
            }
        
        # Compare
        chamber_sum = chamber_result['total_float']
        hodges = hodges_result['amplitude_float']
        
        difference = abs(chamber_sum - hodges)
        relative_diff = difference / abs(hodges) if hodges != 0 else float('inf')
        
        print(f"\n" + "="*80)
        print("COMPARISON")
        print("="*80)
        print(f"Chamber Sum:     {chamber_sum:.10e}")
        print(f"Hodges:          {hodges:.10e}")
        print(f"Difference:      {difference:.10e}")
        print(f"Relative Diff:   {relative_diff:.10e}")
        
        # Check if they match (within numerical tolerance)
        matches = relative_diff < 1e-6
        
        if matches:
            print(f"\n✓ MATCH! Chamber sum equals Hodges amplitude (rel diff < 1e-6)")
            conclusion = "VERIFIED: Gravity amplitude = Sum over all ordering chambers"
        else:
            print(f"\n✗ NO MATCH. Relative difference: {relative_diff:.2e}")
            
            # Check ratio
            ratio = chamber_sum / hodges if hodges != 0 else float('inf')
            print(f"Ratio (chamber/hodges): {ratio:.6f}")
            
            if abs(ratio - 1) < 0.01:
                conclusion = "CLOSE: Within 1% but not exact match"
            elif abs(ratio) > 0.1 and abs(ratio) < 10:
                conclusion = f"PARTIAL: Off by factor of ~{ratio:.2f}"
            else:
                conclusion = f"MISMATCH: Large discrepancy (ratio = {ratio:.2e})"
        
        print(f"\nCONCLUSION: {conclusion}")
        
        return {
            'status': 'ok',
            'matches': matches,
            'chamber_sum': chamber_sum,
            'hodges': hodges,
            'difference': difference,
            'relative_diff': relative_diff,
            'ratio': chamber_sum / hodges if hodges != 0 else None,
            'num_chambers': chamber_result['num_chambers'],
            'conclusion': conclusion
        }
    
    def full_analysis(self):
        """
        Run complete chamber summation analysis.
        
        Returns:
            dict with all results
        """
        print("="*80)
        print("FULL CHAMBER SUMMATION ANALYSIS")
        print("="*80)
        print("Testing hypothesis: Gravity amplitude = Sum over all ordering chambers")
        print("="*80)
        
        # Group solutions by chamber
        self.solve_and_group_by_chamber()
        
        # Verify chamber sum equals Hodges
        result = self.verify_chamber_sum_equals_hodges()
        
        # Final summary
        print("\n" + "="*80)
        print("FINAL SUMMARY")
        print("="*80)
        print(f"Number of scattering equation solutions: {len(self.solutions)}")
        print(f"Number of distinct chambers: {result['num_chambers']}")
        print(f"Chamber sum matches Hodges: {result['matches']}")
        print(f"Conclusion: {result['conclusion']}")
        
        return result


def main(seed=42):
    """
    Main execution function.
    
    Args:
        seed: Random seed for kinematics
        
    Returns:
        ChamberSumAnalysis instance with results
    """
    analysis = ChamberSumAnalysis(seed=seed)
    result = analysis.full_analysis()
    return analysis, result


if __name__ == "__main__":
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    analysis, result = main(seed=seed)

