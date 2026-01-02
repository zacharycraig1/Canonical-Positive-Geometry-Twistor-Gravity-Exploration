#!/usr/bin/env sage
# =============================================================================
# Cross-Ratio Coordinate Analysis for 6-Point MHV Gravity
# =============================================================================
# Tests positivity conditions in gauge-invariant cross-ratio coordinates
# rather than the gauge-dependent z4, z5, z6.
#
# Cross-ratios are SL(2,C) invariant and provide a cleaner description
# of the moduli space M_{0,6}.

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/scattering_solver.sage")
load("src/gravity_proof/psi_matrix.sage")

# Import Python modules
from src.kinematics.spinors import SpinorKinematics


class CrossRatioAnalysis:
    """
    Analyzes scattering equation solutions in cross-ratio coordinates.
    
    Cross-ratios are gauge-invariant and may reveal hidden positivity structure.
    """
    
    def __init__(self, kinematics=None, seed=42):
        """
        Initialize cross-ratio analysis.
        
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
    
    def z_to_cross_ratio(self, z4, z5, z6):
        """
        Convert z-coordinates to cross-ratio coordinates.
        
        With gauge fixing (z1, z2, z3) = (0, 1, infinity):
        
        Standard cross-ratio: (z_i - z_j)(z_k - z_l) / ((z_i - z_k)(z_j - z_l))
        
        We define:
            u = z4 / (z4 - 1)        [related to cross-ratio involving z1, z2, z4]
            v = z5 / (z5 - 1)        [related to cross-ratio involving z1, z2, z5]
            w = z6 / (z6 - 1)        [related to cross-ratio involving z1, z2, z6]
        
        Alternative: more symmetric cross-ratios
            u' = (z4 - 0)(1 - inf) / ((z4 - inf)(0 - 1)) = z4  [simplified]
            
        Actually, use proper cross-ratios of the 6 points:
            Let's use (z1, z2; z3, zi) type cross-ratios
        
        Args:
            z4, z5, z6: Worldsheet coordinates
            
        Returns:
            (u, v, w) cross-ratio coordinates
        """
        # Simple transformation: u = z / (z - 1)
        # This maps:
        #   z = 0 -> u = 0
        #   z = 1 -> u = infinity
        #   z = infinity -> u = 1
        #   z in (0, 1) -> u in (-inf, 0)
        #   z in (1, inf) -> u in (1, inf)
        #   z in (-inf, 0) -> u in (0, 1)
        
        def transform(z):
            if z == 1:
                return float('inf')
            return z / (z - 1)
        
        z4_r = self._to_real(z4)
        z5_r = self._to_real(z5)
        z6_r = self._to_real(z6)
        
        u = transform(z4_r)
        v = transform(z5_r)
        w = transform(z6_r)
        
        return (u, v, w)
    
    def z_to_dihedral_cross_ratios(self, z4, z5, z6):
        """
        Compute dihedral cross-ratios commonly used in amplitude literature.
        
        For n=6, the 3 independent cross-ratios are often:
            u1 = s_12 * s_45 / (s_123 * s_345)
            u2 = s_23 * s_56 / (s_234 * s_456)
            u3 = s_34 * s_61 / (s_345 * s_561)
        
        But these depend on Mandelstam invariants, not just z's.
        
        For pure worldsheet cross-ratios, use:
            chi_1 = (z14 * z25 * z36) / (z15 * z26 * z34)
            chi_2 = (z12 * z35 * z46) / (z13 * z24 * z56)
            etc.
        
        With gauge (z1=0, z2=1, z3=inf):
            z14 = -z4, z15 = -z5, z16 = -z6
            z24 = 1-z4, z25 = 1-z5, z26 = 1-z6
            z34 = inf, z35 = inf, z36 = inf
            z45 = z4-z5, z46 = z4-z6, z56 = z5-z6
        
        Returns:
            dict with various cross-ratio combinations
        """
        z4_r = self._to_real(z4)
        z5_r = self._to_real(z5)
        z6_r = self._to_real(z6)
        
        # z_ij values with gauge fixing
        z14 = -z4_r
        z15 = -z5_r
        z16 = -z6_r
        z24 = 1 - z4_r
        z25 = 1 - z5_r
        z26 = 1 - z6_r
        z45 = z4_r - z5_r
        z46 = z4_r - z6_r
        z56 = z5_r - z6_r
        
        cross_ratios = {}
        
        # Cross-ratio (z1, z4; z2, z5) = (z14 * z25) / (z15 * z24)
        try:
            cross_ratios['chi_1425'] = (z14 * z25) / (z15 * z24) if z15 * z24 != 0 else float('inf')
        except:
            cross_ratios['chi_1425'] = None
        
        # Cross-ratio (z1, z4; z2, z6) = (z14 * z26) / (z16 * z24)
        try:
            cross_ratios['chi_1426'] = (z14 * z26) / (z16 * z24) if z16 * z24 != 0 else float('inf')
        except:
            cross_ratios['chi_1426'] = None
        
        # Cross-ratio (z1, z5; z2, z6) = (z15 * z26) / (z16 * z25)
        try:
            cross_ratios['chi_1526'] = (z15 * z26) / (z16 * z25) if z16 * z25 != 0 else float('inf')
        except:
            cross_ratios['chi_1526'] = None
        
        # Cross-ratio (z4, z5; z1, z2) = (z45 * z12) / (z41 * z52) = (z45 * 1) / (z4 * (1-z5))
        try:
            cross_ratios['chi_4512'] = (z45 * 1) / (z4_r * (1 - z5_r)) if z4_r * (1 - z5_r) != 0 else float('inf')
        except:
            cross_ratios['chi_4512'] = None
        
        # Cross-ratio (z4, z6; z1, z2) = (z46 * 1) / (z4 * (1-z6))
        try:
            cross_ratios['chi_4612'] = (z46 * 1) / (z4_r * (1 - z6_r)) if z4_r * (1 - z6_r) != 0 else float('inf')
        except:
            cross_ratios['chi_4612'] = None
        
        # Cross-ratio (z5, z6; z1, z2) = (z56 * 1) / (z5 * (1-z6))
        try:
            cross_ratios['chi_5612'] = (z56 * 1) / (z5_r * (1 - z6_r)) if z5_r * (1 - z6_r) != 0 else float('inf')
        except:
            cross_ratios['chi_5612'] = None
        
        # Simple ratios u, v, w
        try:
            cross_ratios['u'] = z4_r / (z4_r - 1) if z4_r != 1 else float('inf')
            cross_ratios['v'] = z5_r / (z5_r - 1) if z5_r != 1 else float('inf')
            cross_ratios['w'] = z6_r / (z6_r - 1) if z6_r != 1 else float('inf')
        except:
            pass
        
        return cross_ratios
    
    def _to_real(self, z):
        """Convert to real float."""
        if hasattr(z, 'real') and callable(z.real):
            return float(z.real())
        elif hasattr(z, 'real'):
            return float(z.real)
        else:
            return float(z)
    
    def analyze_solutions_in_cross_ratios(self):
        """
        Analyze all scattering equation solutions in cross-ratio coordinates.
        
        Returns:
            dict with cross-ratio analysis for each solution
        """
        print("="*80)
        print("CROSS-RATIO ANALYSIS")
        print("="*80)
        
        # Solve scattering equations
        print("\nSolving scattering equations...")
        self.solutions = self.solver.solve_numerical()
        print(f"Found {len(self.solutions)} solutions")
        
        results = []
        
        print("\n" + "="*80)
        print("SOLUTIONS IN CROSS-RATIO COORDINATES")
        print("="*80)
        
        for i, sol in enumerate(self.solutions):
            z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
            z4_r = self._to_real(z4)
            z5_r = self._to_real(z5)
            z6_r = self._to_real(z6)
            
            # Simple cross-ratios
            u, v, w = self.z_to_cross_ratio(z4, z5, z6)
            
            # Extended cross-ratios
            cross_ratios = self.z_to_dihedral_cross_ratios(z4, z5, z6)
            
            # Check positivity conditions
            u_positive = u > 0 if u != float('inf') else False
            v_positive = v > 0 if v != float('inf') else False
            w_positive = w > 0 if w != float('inf') else False
            all_positive = u_positive and v_positive and w_positive
            
            u_in_01 = 0 < u < 1 if u != float('inf') else False
            v_in_01 = 0 < v < 1 if v != float('inf') else False
            w_in_01 = 0 < w < 1 if w != float('inf') else False
            all_in_01 = u_in_01 and v_in_01 and w_in_01
            
            result = {
                'index': i,
                'z4': z4_r,
                'z5': z5_r,
                'z6': z6_r,
                'u': u,
                'v': v,
                'w': w,
                'cross_ratios': cross_ratios,
                'u_positive': u_positive,
                'v_positive': v_positive,
                'w_positive': w_positive,
                'all_uvw_positive': all_positive,
                'all_uvw_in_01': all_in_01
            }
            results.append(result)
            
            print(f"\nSolution {i+1}:")
            print(f"  z-coordinates: z4={z4_r:+.6f}, z5={z5_r:+.6f}, z6={z6_r:+.6f}")
            print(f"  Cross-ratios:  u={u:+.6f}, v={v:+.6f}, w={w:+.6f}")
            print(f"  u > 0: {u_positive}, v > 0: {v_positive}, w > 0: {w_positive}")
            print(f"  All (u,v,w) > 0: {all_positive}")
            print(f"  All (u,v,w) in (0,1): {all_in_01}")
        
        return results
    
    def check_positivity_in_cross_ratios(self):
        """
        Check if any positivity condition holds in cross-ratio coordinates.
        
        Returns:
            dict with positivity check results
        """
        results = self.analyze_solutions_in_cross_ratios()
        
        print("\n" + "="*80)
        print("POSITIVITY CHECK SUMMARY")
        print("="*80)
        
        # Count how many solutions satisfy each condition
        count_all_positive = sum(1 for r in results if r['all_uvw_positive'])
        count_all_in_01 = sum(1 for r in results if r['all_uvw_in_01'])
        
        print(f"\nSolutions with all (u,v,w) > 0: {count_all_positive}/{len(results)}")
        print(f"Solutions with all (u,v,w) in (0,1): {count_all_in_01}/{len(results)}")
        
        # Check ordering conditions
        orderings = {
            'u > v > w': 0,
            'u > w > v': 0,
            'v > u > w': 0,
            'v > w > u': 0,
            'w > u > v': 0,
            'w > v > u': 0
        }
        
        for r in results:
            u, v, w = r['u'], r['v'], r['w']
            if u > v > w:
                orderings['u > v > w'] += 1
            elif u > w > v:
                orderings['u > w > v'] += 1
            elif v > u > w:
                orderings['v > u > w'] += 1
            elif v > w > u:
                orderings['v > w > u'] += 1
            elif w > u > v:
                orderings['w > u > v'] += 1
            elif w > v > u:
                orderings['w > v > u'] += 1
        
        print("\nOrderings in (u,v,w) space:")
        for ordering, count in orderings.items():
            if count > 0:
                print(f"  {ordering}: {count} solution(s)")
        
        # Check if any single ordering contains all solutions
        single_ordering = None
        for ordering, count in orderings.items():
            if count == len(results):
                single_ordering = ordering
                break
        
        print("\n" + "="*80)
        print("CONCLUSION")
        print("="*80)
        
        if single_ordering:
            print(f"\n✓ All {len(results)} solutions satisfy ordering: {single_ordering}")
            print("  This could define the positive region in cross-ratio coordinates!")
            conclusion = f"POSITIVE REGION FOUND: {single_ordering}"
        elif count_all_positive == len(results):
            print(f"\n✓ All {len(results)} solutions have u,v,w > 0")
            print("  Positivity is satisfied but no single ordering")
            conclusion = "PARTIAL: All positive but different orderings"
        elif count_all_in_01 == len(results):
            print(f"\n✓ All {len(results)} solutions have u,v,w in (0,1)")
            conclusion = "PARTIAL: All in unit interval"
        else:
            print(f"\n✗ No universal positivity condition found")
            print(f"  Only {count_all_positive}/{len(results)} have all positive cross-ratios")
            conclusion = "NO UNIVERSAL POSITIVITY IN CROSS-RATIOS"
        
        return {
            'results': results,
            'count_all_positive': count_all_positive,
            'count_all_in_01': count_all_in_01,
            'orderings': orderings,
            'single_ordering': single_ordering,
            'conclusion': conclusion
        }
    
    def full_analysis(self):
        """
        Run complete cross-ratio analysis.
        
        Returns:
            dict with all results
        """
        return self.check_positivity_in_cross_ratios()


def main(seed=42):
    """
    Main execution function.
    
    Args:
        seed: Random seed for kinematics
        
    Returns:
        CrossRatioAnalysis instance with results
    """
    analysis = CrossRatioAnalysis(seed=seed)
    result = analysis.full_analysis()
    return analysis, result


if __name__ == "__main__":
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    analysis, result = main(seed=seed)

