#!/usr/bin/env sage
# =============================================================================
# Alternative Gauge Analysis for 6-Point MHV Gravity
# =============================================================================
# Tests whether different gauge choices for (z_a, z_b, z_c) = (0, 1, infinity)
# might place all scattering equation solutions in a positive region.
#
# The standard gauge fixes z1=0, z2=1, z3=infinity. But we can try:
# - Different orderings of (0, 1, infinity)
# - Different choices of which particles to fix

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Import Python modules
from src.kinematics.spinors import SpinorKinematics


class GaugeAnalysis:
    """
    Analyzes scattering equation solutions under different gauge choices.
    """
    
    def __init__(self, kinematics=None, seed=42):
        """
        Initialize gauge analysis.
        
        Args:
            kinematics: Optional SpinorKinematics object
            seed: Random seed for kinematics generation
        """
        if kinematics is None:
            self.kin = SpinorKinematics.random_rational(6, seed=seed)
        else:
            self.kin = kinematics
    
    def get_mandelstam(self, i, j):
        """Get Mandelstam invariant s_ij."""
        return self.kin.angle(i, j) * self.kin.square(i, j)
    
    def scattering_equations_general_gauge(self, fixed_particles, fixed_values, free_vars):
        """
        Construct scattering equations with general gauge fixing.
        
        Args:
            fixed_particles: tuple of 3 particle indices (0-indexed)
            fixed_values: tuple of 3 values (e.g., (0, 1, infinity))
            free_vars: list of 3 Sage variables for the free particles
            
        Returns:
            List of 3 polynomial equations
        """
        # Build z_i values
        z = [None] * 6
        
        # Set fixed values
        for p, v in zip(fixed_particles, fixed_values):
            z[p] = v
        
        # Set free variables
        free_particles = [i for i in range(6) if i not in fixed_particles]
        for i, p in enumerate(free_particles):
            z[p] = free_vars[i]
        
        # Build scattering equations for free particles
        equations = []
        
        for p in free_particles:
            # f_p = sum_{q != p} s_pq / (z_p - z_q) = 0
            eq = 0
            z_p = z[p]
            
            for q in range(6):
                if q == p:
                    continue
                z_q = z[q]
                s_pq = self.get_mandelstam(p, q)
                
                if z_q == oo:
                    # Term vanishes when z_q = infinity
                    continue
                
                eq += s_pq / (z_p - z_q)
            
            equations.append(eq)
        
        return equations, z, free_particles
    
    def solve_with_gauge(self, fixed_particles, fixed_values):
        """
        Solve scattering equations with specified gauge fixing.
        
        Args:
            fixed_particles: tuple of 3 particle indices (0-indexed)
            fixed_values: tuple of 3 values
            
        Returns:
            List of solution dicts
        """
        # Create free variables
        free_particles = [i for i in range(6) if i not in fixed_particles]
        var_names = [f'z{p+1}' for p in free_particles]  # Use 1-indexed names
        free_vars = [var(name) for name in var_names]
        
        # Get equations
        equations, z_all, _ = self.scattering_equations_general_gauge(
            fixed_particles, fixed_values, free_vars
        )
        
        # Clear denominators and solve
        solutions = []
        
        try:
            # Convert to polynomial system
            R = PolynomialRing(QQ, var_names)
            polys = []
            
            for eq in equations:
                # Multiply by all denominators to clear
                eq_simplified = eq.simplify_full()
                # Extract numerator
                numer = eq_simplified.numerator()
                polys.append(R(numer))
            
            # Solve using variety
            I = R.ideal(polys)
            
            # Try to get solutions
            try:
                variety = I.variety(QQbar)
                
                for sol_dict in variety:
                    sol = {}
                    for i, p in enumerate(free_particles):
                        var_name = f'z{p+1}'
                        sol[var_name] = complex(sol_dict[R.gen(i)])
                    solutions.append(sol)
                    
            except Exception as e:
                # Try numerical approach
                print(f"    Variety computation failed: {e}")
                # Fall back to numerical solving
                pass
                
        except Exception as e:
            print(f"    Equation setup failed: {e}")
        
        return solutions, free_particles
    
    def check_positivity(self, solutions, free_particles):
        """
        Check if all solutions have positive free coordinates.
        
        Args:
            solutions: List of solution dicts
            free_particles: List of free particle indices
            
        Returns:
            dict with positivity results
        """
        all_positive = True
        all_ordered = True
        details = []
        
        for sol in solutions:
            vals = []
            for p in free_particles:
                var_name = f'z{p+1}'
                val = sol.get(var_name, None)
                if val is not None:
                    # Get real part
                    if hasattr(val, 'real'):
                        val_r = float(val.real)
                    else:
                        val_r = float(val)
                    vals.append(val_r)
                else:
                    vals.append(None)
            
            # Check positivity
            is_positive = all(v is not None and v > 0 for v in vals)
            
            # Check ordering (descending)
            is_ordered = False
            if all(v is not None for v in vals):
                is_ordered = vals[0] > vals[1] > vals[2] > 0
            
            if not is_positive:
                all_positive = False
            if not is_ordered:
                all_ordered = False
            
            details.append({
                'values': vals,
                'is_positive': is_positive,
                'is_ordered': is_ordered
            })
        
        return {
            'all_positive': all_positive,
            'all_ordered': all_ordered,
            'details': details,
            'num_solutions': len(solutions)
        }
    
    def test_gauge_permutations(self):
        """
        Test all permutations of which particles are fixed and to what values.
        
        For 6 particles, there are C(6,3) = 20 ways to choose which 3 to fix.
        For each choice, there are 3! = 6 ways to assign (0, 1, infinity).
        Total: 120 gauge choices.
        
        We test a subset of interesting ones.
        
        Returns:
            dict with results for each gauge
        """
        print("="*80)
        print("GAUGE ANALYSIS")
        print("="*80)
        print("Testing different gauge choices to find positivity")
        print("="*80)
        
        from itertools import combinations, permutations
        
        results = {}
        
        # Test primary gauges: fix particles (0,1,2), (0,1,5), (0,2,5), etc.
        # Values always (0, 1, infinity)
        
        gauge_choices = [
            # Standard gauge
            ((0, 1, 2), (0, 1, oo), "Standard: z1=0, z2=1, z3=inf"),
            
            # Fix different particles
            ((0, 1, 5), (0, 1, oo), "z1=0, z2=1, z6=inf"),
            ((0, 2, 5), (0, 1, oo), "z1=0, z3=1, z6=inf"),
            ((1, 2, 5), (0, 1, oo), "z2=0, z3=1, z6=inf"),
            
            # Reorder fixed values
            ((0, 1, 2), (oo, 0, 1), "z1=inf, z2=0, z3=1"),
            ((0, 1, 2), (1, oo, 0), "z1=1, z2=inf, z3=0"),
            
            # Fix particles with different helicities
            ((0, 3, 4), (0, 1, oo), "z1=0, z4=1, z5=inf (mixed helicity)"),
            ((2, 3, 4), (0, 1, oo), "z3=0, z4=1, z5=inf (positive helicity)"),
        ]
        
        for fixed_particles, fixed_values, description in gauge_choices:
            print(f"\n{'='*60}")
            print(f"Gauge: {description}")
            print(f"Fixed: particles {tuple(p+1 for p in fixed_particles)} = {fixed_values}")
            print(f"{'='*60}")
            
            try:
                solutions, free_particles = self.solve_with_gauge(fixed_particles, fixed_values)
                
                if len(solutions) == 0:
                    print(f"  No solutions found (may need numerical solver)")
                    results[description] = {'status': 'no_solutions'}
                    continue
                
                print(f"  Found {len(solutions)} solutions")
                
                # Check positivity
                pos_result = self.check_positivity(solutions, free_particles)
                
                print(f"  All positive: {pos_result['all_positive']}")
                print(f"  All ordered: {pos_result['all_ordered']}")
                
                # Show solution details
                for i, detail in enumerate(pos_result['details']):
                    vals = detail['values']
                    val_str = ', '.join([f'{v:.4f}' if v is not None else 'None' for v in vals])
                    pos_str = '✓' if detail['is_positive'] else '✗'
                    ord_str = '✓' if detail['is_ordered'] else '✗'
                    print(f"    Sol {i+1}: [{val_str}] pos:{pos_str} ord:{ord_str}")
                
                results[description] = {
                    'status': 'ok',
                    'num_solutions': len(solutions),
                    'all_positive': pos_result['all_positive'],
                    'all_ordered': pos_result['all_ordered'],
                    'fixed_particles': fixed_particles,
                    'fixed_values': fixed_values
                }
                
                if pos_result['all_positive']:
                    print(f"\n  ✓ ALL SOLUTIONS POSITIVE!")
                if pos_result['all_ordered']:
                    print(f"\n  ✓ ALL SOLUTIONS ORDERED!")
                    
            except Exception as e:
                print(f"  Error: {e}")
                results[description] = {'status': 'error', 'error': str(e)}
        
        # Summary
        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        
        positive_gauges = [k for k, v in results.items() if v.get('all_positive', False)]
        ordered_gauges = [k for k, v in results.items() if v.get('all_ordered', False)]
        
        if positive_gauges:
            print(f"\nGauges with all positive solutions:")
            for g in positive_gauges:
                print(f"  ✓ {g}")
        else:
            print(f"\nNo gauge found with all positive solutions")
        
        if ordered_gauges:
            print(f"\nGauges with all ordered solutions:")
            for g in ordered_gauges:
                print(f"  ✓ {g}")
        else:
            print(f"\nNo gauge found with all ordered solutions")
        
        return results
    
    def full_analysis(self):
        """
        Run complete gauge analysis.
        
        Returns:
            dict with all results
        """
        return self.test_gauge_permutations()


def main(seed=42):
    """
    Main execution function.
    
    Args:
        seed: Random seed for kinematics
        
    Returns:
        GaugeAnalysis instance with results
    """
    analysis = GaugeAnalysis(seed=seed)
    result = analysis.full_analysis()
    return analysis, result


if __name__ == "__main__":
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    analysis, result = main(seed=seed)

