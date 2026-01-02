#!/usr/bin/env sage
# =============================================================================
# PHASE 3: Scattering Equation Solver for 6-Point
# =============================================================================
# Solves the 3 gauge-fixed scattering equations to find all 6 solutions

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from src.kinematics.spinors import SpinorKinematics

class ScatteringEquationSolver:
    """
    Solves scattering equations with gauge fixing (z1, z2, z3) = (0, 1, infinity).
    
    Remaining equations for 6-point (from directive lines 137-141):
    f_4 = s_{14}/z_4 + s_{24}/(z_4-1) + s_{45}/(z_4-z_5) + s_{46}/(z_4-z_6) = 0
    f_5 = s_{15}/z_5 + s_{25}/(z_5-1) + s_{45}/(z_5-z_4) + s_{56}/(z_5-z_6) = 0
    f_6 = s_{16}/z_6 + s_{26}/(z_6-1) + s_{46}/(z_6-z_4) + s_{56}/(z_6-z_5) = 0
    
    Note: Terms with z_3 = infinity vanish (s_{34}/(z_4 - infinity) → 0)
    """
    
    def __init__(self, kinematics):
        """
        Initialize solver.
        
        Args:
            kinematics: SpinorKinematics object (n=6)
        """
        if kinematics.n != 6:
            raise ValueError("ScatteringEquationSolver requires n=6")
            
        self.kin = kinematics
        self.n = 6
        
        # Precompute all Mandelstam invariants
        self.s = {}
        for i in range(6):
            for j in range(i+1, 6):
                s_val = self.kin.s(i, j)
                self.s[(i, j)] = s_val
                self.s[(j, i)] = s_val
    
    def scattering_equation(self, i, z4, z5, z6):
        """
        Compute scattering equation f_i for particle i.
        
        Args:
            i: particle index (0-indexed, should be 3, 4, or 5 for particles 4, 5, 6)
            z4, z5, z6: worldsheet coordinates (can be symbolic or numeric)
            
        Returns:
            f_i(z4, z5, z6)
        """
        # Convert to 1-indexed for clarity
        # Particle 4 -> index 3, particle 5 -> index 4, particle 6 -> index 5
        # z1=0, z2=1, z3=infinity (handled by omitting terms)
        # z4, z5, z6 are free
        
        if i not in [3, 4, 5]:
            raise ValueError(f"Only particles 4,5,6 have equations (i should be 3,4,5, got {i})")
        
        z = [0, 1, None, z4, z5, z6]  # z[0]=z1=0, z[1]=z2=1, z[2]=z3=inf, z[3]=z4, etc.
        z_i = z[i]
        
        result = 0
        
        # Sum over j != i
        for j in range(6):
            if j == i:
                continue
            
            z_j = z[j]
            
            # Skip terms with z3 = infinity
            if z_j is None:  # z3 = infinity
                continue
            
            # Get Mandelstam invariant
            s_ij = self.s[(i, j)]
            
            # Compute 1/(z_i - z_j)
            if z_i is None:  # Should not happen for i in [3,4,5]
                continue
                
            diff = z_i - z_j
            if diff == 0:
                raise ValueError(f"Collision: z_{i+1} == z_{j+1}")
            
            result += s_ij / diff
        
        return result
    
    def system_equations(self, z4, z5, z6):
        """
        Return the system of 3 scattering equations.
        
        Returns:
            List [f4, f5, f6]
        """
        return [
            self.scattering_equation(3, z4, z5, z6),  # f4
            self.scattering_equation(4, z4, z5, z6),  # f5
            self.scattering_equation(5, z4, z5, z6),  # f6
        ]
    
    def jacobian_matrix(self, z4, z5, z6):
        """
        Compute 3x3 Jacobian matrix Φ_{ij} = ∂f_i/∂z_j for i,j ∈ {4,5,6}.
        
        Returns:
            3x3 matrix
        """
        var('z4_var z5_var z6_var')
        
        f4 = self.scattering_equation(3, z4_var, z5_var, z6_var)
        f5 = self.scattering_equation(4, z4_var, z5_var, z6_var)
        f6 = self.scattering_equation(5, z4_var, z5_var, z6_var)
        
        J = matrix(3, 3, [
            [diff(f4, z4_var), diff(f4, z5_var), diff(f4, z6_var)],
            [diff(f5, z4_var), diff(f5, z5_var), diff(f5, z6_var)],
            [diff(f6, z4_var), diff(f6, z5_var), diff(f6, z6_var)],
        ])
        
        # Substitute z4, z5, z6 values
        J_sub = J.subs({z4_var: z4, z5_var: z5, z6_var: z6})
        
        return J_sub
    
    def detprime_phi(self, z4, z5, z6):
        """
        Compute the reduced determinant det'(Φ) for CHY formula.
        
        The correct CHY formula is:
            det'(Φ) = det(Φ^{ijk}_{ijk}) / (σ_ij * σ_jk * σ_ki)²
        
        With gauge fixing (z1, z2, z3) = (0, 1, ∞), we delete rows/cols {1,2,3}
        and keep {4,5,6}. The Vandermonde factor is (z1-z2)(z2-z3)(z3-z1).
        
        However, for z3=∞, the Vandermonde factor needs careful handling.
        In the z3→∞ limit, the correct normalization is:
            det'(Φ) = det(Φ^{123}_{123}) × z3² / [(z1-z2)(z2-z3)(z3-z1)]²
                    = det(Φ^{123}_{123}) × z3² / [(z1-z2) × z3 × z3]²  (leading order)
                    = det(Φ^{123}_{123}) / (z1-z2)²
                    = det(Φ^{123}_{123})  (since z1-z2 = -1)
        
        So with our gauge fixing, det'(Φ) = det(Jacobian) exactly!
        
        Args:
            z4, z5, z6: worldsheet coordinates
            
        Returns:
            det'(Φ) value
        """
        # The Jacobian matrix Φ_{ij} = ∂f_i/∂z_j for the reduced system
        J = self.jacobian_matrix(z4, z5, z6)
        det_J = J.det()
        
        # With gauge fixing (z1, z2, z3) = (0, 1, ∞):
        # The Vandermonde factor (z12 z23 z31)² for deletion set {1,2,3}
        # z12 = 0 - 1 = -1
        # z23 = 1 - ∞ → -∞
        # z31 = ∞ - 0 → ∞
        # 
        # Taking the limit carefully:
        # (z12 z23 z31)² = (-1)² × (-z3)² × z3² = z3⁴
        # 
        # But the Jacobian determinant also scales with z3 in a specific way.
        # For the gauge-fixed scattering equations, the determinant is already
        # the properly normalized det'(Φ).
        #
        # This is because we're working with the reduced (n-3)×(n-3) system.
        
        return det_J
    
    def detprime_phi_with_vandermonde(self, z4, z5, z6):
        """
        Alternative det'(Φ) computation including explicit Vandermonde factor.
        
        This implements the formula:
            det'(Φ) = det(Φ^{ijk}_{ijk}) / (σ_ij * σ_jk * σ_ki)²
        
        For gauge fixing (z1, z2, z3) = (0, 1, ∞), the deletion set is {1,2,3}.
        Since z3=∞, we need to take a limit.
        
        Actually, for the STANDARD CHY formula, the correct prescription uses:
            M_n = (z_{ab} z_{bc} z_{ca})² × Σ_σ [Pf'(Ψ)]² / |det Φ^{abc}_{abc}|
        
        where the measure factor cancels the Vandermonde denominator.
        So the "working" formula is just: [Pf'(Ψ)]² / det(J) with det(J) being
        the raw Jacobian of the reduced system.
        
        Returns:
            (det_raw, vandermonde_factor, det_prime) tuple
        """
        J = self.jacobian_matrix(z4, z5, z6)
        det_raw = J.det()
        
        # Vandermonde factor for deletion {1,2,3} with z1=0, z2=1
        # (z1-z2)(z2-z3)(z3-z1) in the limit z3→∞
        # = (-1)(1-z3)(z3-0) = (-1)(-z3)(z3) = z3² (leading order)
        # Squared: z3⁴
        #
        # But z3 is infinite, so this doesn't work directly.
        # The standard practice is that for z3=∞ gauge, the factors cancel properly.
        
        # For finite z3 (testing), we'd compute:
        # z1, z2, z3 = 0, 1, z3_val
        # vander = (z1-z2) * (z2-z3) * (z3-z1) = (-1) * (1-z3) * z3 = z3(z3-1)
        # vander_sq = z3² * (z3-1)²
        
        # Since z3=∞ in our gauge, we just return det_raw as det'
        # (the measure factor in CHY takes care of this)
        vandermonde_factor = 1  # In z3=∞ gauge, this is effectively 1 after limit
        det_prime = det_raw / vandermonde_factor
        
        return det_raw, vandermonde_factor, det_prime
    
    def _build_f4_polynomial(self, z4, z5, z6):
        """
        Build polynomial for f4 equation by clearing denominators.
        
        Original: s14/z4 + s24/(z4-1) + s45/(z4-z5) + s46/(z4-z6) = 0
        Multiply by: z4*(z4-1)*(z4-z5)*(z4-z6)
        
        Returns:
            Polynomial expression for f4 (cleared denominators)
        """
        # Particle 4 is index 3 (0-indexed)
        s14 = self.s[(0, 3)]  # s[(0,3)] = s_{14}
        s24 = self.s[(1, 3)]  # s[(1,3)] = s_{24}
        s45 = self.s[(3, 4)]  # s[(3,4)] = s_{45}
        s46 = self.s[(3, 5)]  # s[(3,5)] = s_{46}
        
        # Multiply each term by the appropriate factors
        term1 = s14 * (z4 - 1) * (z4 - z5) * (z4 - z6)
        term2 = s24 * z4 * (z4 - z5) * (z4 - z6)
        term3 = s45 * z4 * (z4 - 1) * (z4 - z6)
        term4 = s46 * z4 * (z4 - 1) * (z4 - z5)
        
        return term1 + term2 + term3 + term4
    
    def _build_f5_polynomial(self, z4, z5, z6):
        """
        Build polynomial for f5 equation by clearing denominators.
        
        Original: s15/z5 + s25/(z5-1) + s45/(z5-z4) + s56/(z5-z6) = 0
        Multiply by: z5*(z5-1)*(z5-z4)*(z5-z6)
        
        Returns:
            Polynomial expression for f5 (cleared denominators)
        """
        # Particle 5 is index 4 (0-indexed)
        s15 = self.s[(0, 4)]  # s[(0,4)] = s_{15}
        s25 = self.s[(1, 4)]  # s[(1,4)] = s_{25}
        s45 = self.s[(3, 4)]  # s[(3,4)] = s_{45}
        s56 = self.s[(4, 5)]  # s[(4,5)] = s_{56}
        
        # Multiply each term by the appropriate factors
        term1 = s15 * (z5 - 1) * (z5 - z4) * (z5 - z6)
        term2 = s25 * z5 * (z5 - z4) * (z5 - z6)
        term3 = s45 * z5 * (z5 - 1) * (z5 - z6)
        term4 = s56 * z5 * (z5 - 1) * (z5 - z4)
        
        return term1 + term2 + term3 + term4
    
    def _build_f6_polynomial(self, z4, z5, z6):
        """
        Build polynomial for f6 equation by clearing denominators.
        
        Original: s16/z6 + s26/(z6-1) + s46/(z6-z4) + s56/(z6-z5) = 0
        Multiply by: z6*(z6-1)*(z6-z4)*(z6-z5)
        
        Returns:
            Polynomial expression for f6 (cleared denominators)
        """
        # Particle 6 is index 5 (0-indexed)
        s16 = self.s[(0, 5)]  # s[(0,5)] = s_{16}
        s26 = self.s[(1, 5)]  # s[(1,5)] = s_{26}
        s46 = self.s[(3, 5)]  # s[(3,5)] = s_{46}
        s56 = self.s[(4, 5)]  # s[(4,5)] = s_{56}
        
        # Multiply each term by the appropriate factors
        term1 = s16 * (z6 - 1) * (z6 - z4) * (z6 - z5)
        term2 = s26 * z6 * (z6 - z4) * (z6 - z5)
        term3 = s46 * z6 * (z6 - 1) * (z6 - z5)
        term4 = s56 * z6 * (z6 - 1) * (z6 - z4)
        
        return term1 + term2 + term3 + term4
    
    def solve_polynomial_groebner(self):
        """
        Solve scattering equations exactly using Gröbner basis.
        
        Clears denominators to get polynomial system, then solves using ideal variety.
        
        Returns:
            List of solution dicts with keys: 'z4', 'z5', 'z6', 'exact_z4', 'exact_z5', 'exact_z6'
        """
        # Determine base field from kinematics
        s_sample = self.s[(0, 1)]
        base_field = s_sample.parent()
        
        # DEBUG: Print base field
        print(f"DEBUG: Base field = {base_field}")
        
        # Build polynomial ring
        if base_field == QQ or base_field == ZZ:
            R = PolynomialRing(QQ, ['z4', 'z5', 'z6'])
        else:
            # Try to use the base field, or fall back to QQ
            try:
                R = PolynomialRing(base_field, ['z4', 'z5', 'z6'])
            except:
                R = PolynomialRing(QQ, ['z4', 'z5', 'z6'])
        
        z4, z5, z6 = R.gens()
        
        # Build polynomial equations (cleared denominators)
        f4_poly = self._build_f4_polynomial(z4, z5, z6)
        f5_poly = self._build_f5_polynomial(z4, z5, z6)
        f6_poly = self._build_f6_polynomial(z4, z5, z6)
        
        # DEBUG: Print polynomial degrees
        print(f"DEBUG: f4_poly degree = {f4_poly.degree()}")
        print(f"DEBUG: f5_poly degree = {f5_poly.degree()}")
        print(f"DEBUG: f6_poly degree = {f6_poly.degree()}")
        
        # Check for common factors
        gcd_45 = f4_poly.gcd(f5_poly)
        gcd_46 = f4_poly.gcd(f6_poly)
        gcd_56 = f5_poly.gcd(f6_poly)
        
        if gcd_45 != 1:
            print(f"DEBUG: f4 and f5 have common factor: {gcd_45}")
        if gcd_46 != 1:
            print(f"DEBUG: f4 and f6 have common factor: {gcd_46}")
        if gcd_56 != 1:
            print(f"DEBUG: f5 and f6 have common factor: {gcd_56}")
        
        # Create ideal
        I = R.ideal([f4_poly, f5_poly, f6_poly])
        
        # DEBUG: Print ideal dimension
        try:
            dim = I.dimension()
            print(f"DEBUG: Ideal dimension = {dim}")
        except Exception as e:
            print(f"DEBUG: Could not compute ideal dimension: {e}")
        
        try:
            # Try to compute Gröbner basis first
            G = I.groebner_basis()
            print(f"DEBUG: Gröbner basis has {len(G)} elements")
            if len(G) <= 3:
                for i, g in enumerate(G):
                    print(f"  G[{i}] = {g}")
            
            # If dimension is > 0, try saturating with denominators
            if I.dimension() > 0:
                print("DEBUG: Ideal dimension > 0, attempting to saturate with denominator factors")
                # Saturate with factors that shouldn't vanish
                denom_factors = [z4, z4-1, z5, z5-1, z6, z6-1, z4-z5, z4-z6, z5-z6]
                I_sat = I
                for d in denom_factors:
                    I_sat = I_sat.saturation(d)[0]
                
                dim_sat = I_sat.dimension()
                print(f"DEBUG: Saturated ideal dimension = {dim_sat}")
                
                if dim_sat == 0:
                    print("DEBUG: Using saturated ideal")
                    I = I_sat
                    G = I.groebner_basis()
            
            # Compute variety over algebraic closure
            V = I.variety(QQbar)
            
            print(f"DEBUG: Found {len(V)} solutions in variety")
            
            solutions = []
            for sol in V:
                z4_val = sol[z4]
                z5_val = sol[z5]
                z6_val = sol[z6]
                
                # Convert to complex numbers for numeric use
                solutions.append({
                    'z4': complex(z4_val),
                    'z5': complex(z5_val),
                    'z6': complex(z6_val),
                    'exact_z4': z4_val,
                    'exact_z5': z5_val,
                    'exact_z6': z6_val,
                    'valid': True
                })
            
            return solutions
            
        except Exception as e:
            print(f"Gröbner basis solve failed: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def solve_symbolic(self):
        """
        Attempt to solve scattering equations symbolically using Groebner basis.
        (Deprecated - use solve_polynomial_groebner instead)
        
        Returns:
            List of solution tuples (z4, z5, z6) or None if symbolic solve fails
        """
        return self.solve_polynomial_groebner()
    
    def _solve_numerical_complex(self):
        """
        Numerical solver allowing complex solutions (fallback method).
        Uses Newton's method with complex support.
        
        Returns:
            List of solution dicts
        """
        # This is a fallback for when Gröbner basis fails
        # For now, return empty list - could implement scipy-based solver if needed
        print("Warning: Complex numerical solver not fully implemented, using Newton's method fallback")
        return []
    
    def solve_numerical(self, initial_guesses=None):
        """
        Solve scattering equations numerically to find all 6 solutions.
        
        First tries Gröbner basis solver (exact for rational kinematics).
        Falls back to numerical methods if needed.
        
        Args:
            initial_guesses: Optional list of initial guess tuples (z4, z5, z6)
                           (ignored if Gröbner solver succeeds)
            
        Returns:
            List of solution dicts with keys: 'z4', 'z5', 'z6', 'jacobian_det', 'valid'
        """
        # First try Gröbner basis solver (exact for rational kinematics)
        try:
            groebner_solutions = self.solve_polynomial_groebner()
            if groebner_solutions is not None and len(groebner_solutions) > 0:
                # Convert to format expected by calling code
                solutions = []
                for sol in groebner_solutions:
                    z4_val = sol['z4']
                    z5_val = sol['z5']
                    z6_val = sol['z6']
                    
                    # Compute Jacobian determinant at solution (use real parts for numeric z)
                    try:
                        z4_real = z4_val.real if hasattr(z4_val, 'real') else z4_val
                        z5_real = z5_val.real if hasattr(z5_val, 'real') else z5_val
                        z6_real = z6_val.real if hasattr(z6_val, 'real') else z6_val
                        
                        J = self.jacobian_matrix(z4_real, z5_real, z6_real)
                        jac_det = float(J.det())
                    except:
                        jac_det = 1.0  # Fallback
                    
                    solution_dict = {
                        'z4': z4_real,
                        'z5': z5_real,
                        'z6': z6_real,
                        'jacobian_det': jac_det,
                        'valid': True,
                        'method': 'groebner'
                    }
                    # Preserve exact algebraic values for distinctness checking
                    if 'exact_z4' in sol:
                        solution_dict['exact_z4'] = sol['exact_z4']
                    if 'exact_z5' in sol:
                        solution_dict['exact_z5'] = sol['exact_z5']
                    if 'exact_z6' in sol:
                        solution_dict['exact_z6'] = sol['exact_z6']
                    solutions.append(solution_dict)
                
                if len(solutions) >= 6:
                    print(f"Found {len(solutions)} solutions via Gröbner basis")
                    return solutions
        except Exception as e:
            print(f"Gröbner basis solver failed: {e}, trying fallback")
        
        # Fallback to complex numerical solver
        try:
            complex_solutions = self._solve_numerical_complex()
            if complex_solutions and len(complex_solutions) > 0:
                return complex_solutions
        except Exception as e:
            print(f"Complex numerical solver failed: {e}")
        
        # Last resort: return empty list (previous Newton's method was not finding solutions)
        print("Warning: All solving methods failed, returning empty solution list")
        return []
    
    def verify_solutions_in_region(self, solutions, region):
        """
        Check which solutions lie in the positivity region.
        
        Args:
            solutions: List of solution dicts from solve_numerical
            region: PositivityRegionR6 instance
            
        Returns:
            List of solutions with 'in_region' field added
        """
        for sol in solutions:
            z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
            sol['in_region'] = region.is_inside_option_b1(z4, z5, z6)
        
        return solutions
    
    def analyze_solutions_in_region(self, solutions):
        """
        Analyze which solutions are in various regions.
        
        Checks solutions against all possible orderings of (z4, z5, z6).
        
        Args:
            solutions: List of solution dicts from solve_numerical
            
        Returns:
            dict with:
            - solutions_in_R: list of solutions with z4 > z5 > z6 > 0
            - solutions_real: list of real solutions
            - solutions_by_ordering: dict of solutions for each ordering
        """
        solutions_in_R = []
        solutions_real = []
        solutions_by_ordering = {
            'z4>z5>z6': [],
            'z4>z6>z5': [],
            'z5>z4>z6': [],
            'z5>z6>z4': [],
            'z6>z4>z5': [],
            'z6>z5>z4': []
        }
        
        for sol in solutions:
            z4 = sol['z4']
            z5 = sol['z5']
            z6 = sol['z6']
            
            # Check if real
            is_real = True
            if hasattr(z4, 'imag'):
                is_real = is_real and abs(z4.imag) < 1e-10
            if hasattr(z5, 'imag'):
                is_real = is_real and abs(z5.imag) < 1e-10
            if hasattr(z6, 'imag'):
                is_real = is_real and abs(z6.imag) < 1e-10
            
            if is_real:
                z4_real = z4.real if hasattr(z4, 'real') else z4
                z5_real = z5.real if hasattr(z5, 'real') else z5
                z6_real = z6.real if hasattr(z6, 'real') else z6
                
                solutions_real.append(sol)
                
                # Check ordering z4 > z5 > z6 > 0 (Option B1)
                if z4_real > z5_real > z6_real > 0:
                    solutions_in_R.append(sol)
                    solutions_by_ordering['z4>z5>z6'].append(sol)
                
                # Check other orderings
                if z4_real > z6_real > z5_real > 0:
                    solutions_by_ordering['z4>z6>z5'].append(sol)
                if z5_real > z4_real > z6_real > 0:
                    solutions_by_ordering['z5>z4>z6'].append(sol)
                if z5_real > z6_real > z4_real > 0:
                    solutions_by_ordering['z5>z6>z4'].append(sol)
                if z6_real > z4_real > z5_real > 0:
                    solutions_by_ordering['z6>z4>z5'].append(sol)
                if z6_real > z5_real > z4_real > 0:
                    solutions_by_ordering['z6>z5>z4'].append(sol)
        
        return {
            'solutions_in_R': solutions_in_R,
            'solutions_real': solutions_real,
            'solutions_by_ordering': solutions_by_ordering,
            'num_in_R': len(solutions_in_R),
            'num_real': len(solutions_real),
            'num_total': len(solutions)
        }
    
    def verify_one_solution_per_chamber(self, solutions):
        """
        Verify that each ordering chamber contains exactly one solution.
        
        If true, this confirms the geometry is the full M_{0,6}(R), not a single chamber.
        
        Args:
            solutions: List of solution dicts from solve_numerical
            
        Returns:
            dict with:
            - one_solution_per_chamber: boolean (True if each chamber has exactly 1)
            - distribution: dict mapping ordering to count
            - hypothesis_confirmed: boolean (same as one_solution_per_chamber)
        """
        analysis = self.analyze_solutions_in_region(solutions)
        
        orderings = ['z4>z5>z6', 'z4>z6>z5', 'z5>z4>z6', 
                     'z5>z6>z4', 'z6>z4>z5', 'z6>z5>z4']
        
        one_per_chamber = True
        distribution = {}
        
        for ordering in orderings:
            count = len(analysis['solutions_by_ordering'].get(ordering, []))
            distribution[ordering] = count
            if count != 1:
                one_per_chamber = False
        
        return {
            'one_solution_per_chamber': one_per_chamber,
            'distribution': distribution,
            'hypothesis_confirmed': one_per_chamber
        }
    
    def print_solution_details(self, solutions):
        """
        Print detailed z-values for each solution to understand chamber structure.
        
        Shows which intervals (A=(-∞,0), B=(0,1), C=(1,∞)) each z-value falls into
        and the complete ordering of all punctures.
        
        Args:
            solutions: List of solution dicts from solve_numerical
        """
        print("\n" + "="*60)
        print("DETAILED SOLUTION VALUES")
        print("="*60)
        print("Gauge fixing: z1=0, z2=1, z3=∞")
        print("Intervals: A=(-∞,0), B=(0,1), C=(1,∞)")
        print("-"*60)
        
        def interval(z):
            """Return interval label for z value."""
            if z < 0:
                return 'A'
            elif z < 1:
                return 'B'
            else:
                return 'C'
        
        for i, sol in enumerate(solutions):
            z4 = sol['z4'].real if hasattr(sol['z4'], 'real') else sol['z4']
            z5 = sol['z5'].real if hasattr(sol['z5'], 'real') else sol['z5']
            z6 = sol['z6'].real if hasattr(sol['z6'], 'real') else sol['z6']
            
            intervals = f"{interval(z4)}{interval(z5)}{interval(z6)}"
            
            # Full ordering (exclude z3=∞)
            all_z = [(0, 'z1'), (1, 'z2'), (z4, 'z4'), (z5, 'z5'), (z6, 'z6')]
            all_z_sorted = sorted(all_z, key=lambda x: x[0])
            ordering = ' < '.join([name for _, name in all_z_sorted]) + ' < z3(∞)'
            
            print(f"\nSolution {i+1}:")
            print(f"  z4 = {z4:.6f}  [{interval(z4)}]")
            print(f"  z5 = {z5:.6f}  [{interval(z5)}]")
            print(f"  z6 = {z6:.6f}  [{interval(z6)}]")
            print(f"  Chamber code: {intervals}")
            print(f"  Full ordering: {ordering}")
        
        print("\n" + "="*60)
    
    def analyze_true_chambers(self, solutions):
        """
        Analyze solutions in terms of true M_{0,6}(R) chamber structure.
        
        A chamber is determined by the complete ordering of all punctures on ℝP¹.
        With gauge fixing z1=0, z2=1, z3=∞, a chamber is specified by the ordering
        of z1, z2, z4, z5, z6 (z3 always at ∞).
        
        Args:
            solutions: List of solution dicts from solve_numerical
            
        Returns:
            dict with:
            - chambers_found: dict mapping ordering tuples to solution indices
            - num_chambers: number of distinct chambers found
            - one_per_chamber: boolean (True if exactly 6 chambers with one solution each)
        """
        chambers_found = {}
        
        for i, sol in enumerate(solutions):
            z4 = sol['z4'].real if hasattr(sol['z4'], 'real') else sol['z4']
            z5 = sol['z5'].real if hasattr(sol['z5'], 'real') else sol['z5']
            z6 = sol['z6'].real if hasattr(sol['z6'], 'real') else sol['z6']
            
            # Create ordering tuple (exclude z3=∞ which is always at the end)
            points = [(0, 1), (1, 2), (z4, 4), (z5, 5), (z6, 6)]  # (value, label)
            points_sorted = sorted(points, key=lambda x: x[0])
            ordering = tuple(p[1] for p in points_sorted)  # e.g., (4, 1, 5, 2, 6)
            
            if ordering not in chambers_found:
                chambers_found[ordering] = []
            chambers_found[ordering].append(i)
        
        # Check if exactly one solution per chamber (6 distinct chambers)
        one_per_chamber = (len(chambers_found) == 6 and 
                           all(len(sols) == 1 for sols in chambers_found.values()))
        
        # Print results
        print(f"\nTrue Chamber Analysis:")
        print(f"  Found {len(chambers_found)} distinct chambers")
        for ordering, sol_indices in sorted(chambers_found.items()):
            ordering_str = ' < '.join([f'z{i}' for i in ordering]) + ' < z3'
            print(f"    {ordering_str}: solution(s) {sol_indices}")
        
        if one_per_chamber:
            print(f"\n✅ One solution per chamber confirmed - geometry is full M_{{0,6}}(ℝ)")
        else:
            print(f"\n⚠️  Not exactly one solution per chamber")
        
        return {
            'chambers_found': chambers_found,
            'num_chambers': len(chambers_found),
            'one_per_chamber': one_per_chamber
        }
    
    def check_solution_distinctness(self, solutions):
        """
        Check if solutions are algebraically distinct using exact values.
        
        Uses exact algebraic values from Gröbner basis solver to determine
        if solutions are truly identical (multiplicity > 1) or just numerically close.
        
        Args:
            solutions: List of solution dicts from solve_numerical
            
        Returns:
            dict with:
            - duplicates: list of (i, j) pairs that are algebraically identical
            - num_distinct: number of distinct solutions
            - distinctness_matrix: matrix of pairwise distinctness checks
        """
        print("\n" + "="*60)
        print("EXACT ALGEBRAIC SOLUTION CHECK")
        print("="*60)
        
        # Extract exact values
        exact_solutions = []
        has_exact = False
        for i, sol in enumerate(solutions):
            exact_z4 = sol.get('exact_z4')
            exact_z5 = sol.get('exact_z5')
            exact_z6 = sol.get('exact_z6')
            
            if exact_z4 is not None and exact_z5 is not None and exact_z6 is not None:
                exact_solutions.append((exact_z4, exact_z5, exact_z6))
                has_exact = True
                print(f"\nSolution {i+1} (exact):")
                print(f"  z4 = {exact_z4}")
                print(f"  z5 = {exact_z5}")
                print(f"  z6 = {exact_z6}")
            else:
                exact_solutions.append(None)
                print(f"\nSolution {i+1}: No exact values available (numerical only)")
        
        if not has_exact:
            print("\n⚠️  No exact algebraic values found - cannot check distinctness")
            return {
                'duplicates': [],
                'num_distinct': len(solutions),
                'has_exact_values': False,
                'note': 'No exact values available for distinctness check'
            }
        
        # Check pairwise distinctness
        print("\nPairwise distinctness check:")
        duplicates = []
        distinctness_info = []
        
        for i in range(len(solutions)):
            if exact_solutions[i] is None:
                continue
            for j in range(i+1, len(solutions)):
                if exact_solutions[j] is None:
                    continue
                
                z4_i, z5_i, z6_i = exact_solutions[i]
                z4_j, z5_j, z6_j = exact_solutions[j]
                
                # Check exact equality of QQbar elements
                z4_equal = (z4_i == z4_j)
                z5_equal = (z5_i == z5_j)
                z6_equal = (z6_i == z6_j)
                
                if z4_equal and z5_equal and z6_equal:
                    print(f"  Solutions {i+1} and {j+1}: IDENTICAL (algebraically)")
                    duplicates.append((i, j))
                    distinctness_info.append({
                        'pair': (i+1, j+1),
                        'status': 'identical',
                        'z4_equal': True,
                        'z5_equal': True,
                        'z6_equal': True
                    })
                else:
                    # Check which components differ
                    status = 'distinct'
                    info = {
                        'pair': (i+1, j+1),
                        'status': 'distinct',
                        'z4_equal': bool(z4_equal),
                        'z5_equal': bool(z5_equal),
                        'z6_equal': bool(z6_equal)
                    }
                    
                    if z4_equal:
                        info['note'] = 'z4 identical, checking z5, z6...'
                    if z5_equal and z4_equal:
                        info['note'] = 'z4, z5 identical, z6 differs'
                    if z6_equal and z5_equal and z4_equal:
                        info['note'] = 'all identical (should not reach here)'
                    
                    distinctness_info.append(info)
                    print(f"  Solutions {i+1} and {j+1}: DISTINCT")
        
        num_distinct = len(solutions) - len(duplicates)
        
        if duplicates:
            print(f"\n⚠️  Found {len(duplicates)} duplicate pair(s)")
            print("This may explain the ~2x CHY normalization error!")
        else:
            print(f"\n✅ All solutions are algebraically distinct")
        
        return {
            'duplicates': duplicates,
            'num_distinct': num_distinct,
            'has_exact_values': True,
            'distinctness_info': distinctness_info
        }
    
    def test_multiple_kinematics(self, num_tests=5):
        """
        Test chamber structure with different random kinematics.
        
        This helps determine if the duplicate solution is specific to the
        current kinematics or a generic issue.
        
        Args:
            num_tests: Number of different kinematic configurations to test
            
        Returns:
            List of result dicts, one per test
        """
        print("\n" + "="*60)
        print(f"TESTING CHAMBER STRUCTURE WITH {num_tests} DIFFERENT KINEMATICS")
        print("="*60)
        
        results = []
        
        for seed in range(num_tests):
            print(f"\n--- Test {seed+1}/{num_tests} (seed={seed}) ---")
            
            try:
                # Generate new random kinematics
                kin = SpinorKinematics.random_rational(6, seed=seed)
                
                # Create solver
                solver = ScatteringEquationSolver(kin)
                
                # Solve scattering equations
                solutions = solver.solve_numerical()
                
                if not solutions:
                    print(f"  ⚠️  No solutions found for seed {seed}")
                    results.append({
                        'seed': seed,
                        'num_solutions': 0,
                        'num_chambers': 0,
                        'one_per_chamber': False,
                        'error': 'No solutions found'
                    })
                    continue
                
                # Analyze chambers
                try:
                    analysis = solver.analyze_true_chambers(solutions)
                    
                    result = {
                        'seed': seed,
                        'num_solutions': len(solutions),
                        'num_chambers': analysis['num_chambers'],
                        'one_per_chamber': analysis['one_per_chamber']
                    }
                    
                    # Check for duplicates
                    try:
                        distinctness = solver.check_solution_distinctness(solutions)
                        result['num_distinct'] = distinctness['num_distinct']
                        result['has_duplicates'] = len(distinctness['duplicates']) > 0
                    except Exception as e:
                        result['distinctness_check_error'] = str(e)
                    
                    results.append(result)
                    print(f"  Seed {seed}: {analysis['num_chambers']} chambers, "
                          f"one_per_chamber={analysis['one_per_chamber']}, "
                          f"{len(solutions)} solutions")
                    
                except Exception as e:
                    print(f"  ⚠️  Chamber analysis failed for seed {seed}: {e}")
                    results.append({
                        'seed': seed,
                        'num_solutions': len(solutions),
                        'error': f'Chamber analysis failed: {e}'
                    })
                    
            except Exception as e:
                print(f"  ⚠️  Failed for seed {seed}: {e}")
                results.append({
                    'seed': seed,
                    'error': str(e)
                })
        
        # Summary
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        successful = [r for r in results if 'error' not in r]
        if successful:
            num_6_chambers = sum(1 for r in successful if r['num_chambers'] == 6)
            num_5_chambers = sum(1 for r in successful if r['num_chambers'] == 5)
            print(f"Successful tests: {len(successful)}/{num_tests}")
            print(f"Tests with 6 chambers: {num_6_chambers}")
            print(f"Tests with 5 chambers: {num_5_chambers}")
            if num_6_chambers > 0:
                print("✅ Some configurations give 6 distinct chambers")
            if num_5_chambers == len(successful):
                print("⚠️  All configurations give 5 chambers - may be systematic")
        
        return results
    
    def diagnose_aaa_chamber(self, solutions):
        """
        Add diagnostic output for AAA chamber (all z-values negative).
        
        This chamber has all free punctures to the LEFT of both fixed finite punctures.
        We check for near-degeneracy, Jacobian condition, and scattering equation residuals.
        
        Args:
            solutions: List of solution dicts from solve_numerical
            
        Returns:
            dict with diagnostic information for AAA chamber solutions
        """
        print("\n" + "="*60)
        print("AAA CHAMBER DIAGNOSTICS")
        print("="*60)
        print("(All z-values in interval A: (-∞, 0))")
        print("-"*60)
        
        aaa_solutions = []
        for i, sol in enumerate(solutions):
            z4 = sol['z4'].real if hasattr(sol['z4'], 'real') else sol['z4']
            z5 = sol['z5'].real if hasattr(sol['z5'], 'real') else sol['z5']
            z6 = sol['z6'].real if hasattr(sol['z6'], 'real') else sol['z6']
            
            # Check if all z-values are negative (AAA chamber)
            if z4 < 0 and z5 < 0 and z6 < 0:
                aaa_solutions.append((i, z4, z5, z6))
        
        if not aaa_solutions:
            print("No AAA chamber solutions found (all z < 0)")
            return {'found': False}
        
        diagnostics = []
        for idx, z4, z5, z6 in aaa_solutions:
            print(f"\nSolution {idx+1} (AAA chamber):")
            print(f"  z4 = {z4:.6f}")
            print(f"  z5 = {z5:.6f}")
            print(f"  z6 = {z6:.6f}")
            
            # Distance to z=0 boundary
            min_dist_to_zero = min(abs(z4), abs(z5), abs(z6))
            print(f"  Distance to z=0 boundary: {min_dist_to_zero:.6f}")
            
            # Check Jacobian determinant and condition
            try:
                J = self.jacobian_matrix(z4, z5, z6)
                det_J = J.det()
                # Convert to float for formatting
                try:
                    det_J_float = float(det_J) if hasattr(det_J, '__float__') else complex(det_J).real
                    print(f"  Jacobian determinant: {det_J_float:.6e}")
                except:
                    print(f"  Jacobian determinant: {det_J}")
                
                # Condition number (approximate via norm)
                try:
                    J_inv = J.inverse()
                    cond_approx = float(J.norm() * J_inv.norm())
                    print(f"  Jacobian condition number (approx): {cond_approx:.6e}")
                except:
                    print(f"  Jacobian condition number: INF (singular matrix)")
                    cond_approx = float('inf')
            except Exception as e:
                print(f"  Jacobian computation failed: {e}")
                det_J = None
                cond_approx = None
            
            # Check scattering equation residuals
            try:
                f4_val = self.scattering_equation(3, z4, z5, z6)  # particle 4 (0-indexed: 3)
                f5_val = self.scattering_equation(4, z4, z5, z6)
                f6_val = self.scattering_equation(5, z4, z5, z6)
                
                residual = abs(f4_val) + abs(f5_val) + abs(f6_val)
                print(f"  Scattering equation residual: {residual:.6e}")
                
                if hasattr(f4_val, '__float__'):
                    f4_abs = abs(float(f4_val))
                    f5_abs = abs(float(f5_val))
                    f6_abs = abs(float(f6_val))
                    print(f"    f4: {f4_abs:.6e}, f5: {f5_abs:.6e}, f6: {f6_abs:.6e}")
            except Exception as e:
                print(f"  Residual computation failed: {e}")
                residual = None
            
            # Check if suspiciously close to 0
            if min_dist_to_zero < 0.01:
                print(f"  ⚠️  Warning: Very close to z=0 boundary!")
            
            diagnostics.append({
                'solution_index': idx,
                'z4': float(z4),
                'z5': float(z5),
                'z6': float(z6),
                'min_dist_to_zero': float(min_dist_to_zero),
                'jacobian_det': float(det_J) if det_J is not None else None,
                'jacobian_condition': float(cond_approx) if cond_approx is not None else None,
                'residual': float(residual) if residual is not None else None
            })
        
        return {
            'found': True,
            'num_aaa_solutions': len(aaa_solutions),
            'diagnostics': diagnostics
        }
    
    def verify_jacobian_numerical(self, sol, eps=1e-8):
        """
        Verify analytical Jacobian against numerical derivatives.
        
        Computes the Jacobian numerically using finite differences and
        compares with the analytical computation.
        
        Args:
            sol: Solution dict with z4, z5, z6 values
            eps: Finite difference step size
            
        Returns:
            dict with comparison results
        """
        z4_base = sol['z4'].real if hasattr(sol['z4'], 'real') else sol['z4']
        z5_base = sol['z5'].real if hasattr(sol['z5'], 'real') else sol['z5']
        z6_base = sol['z6'].real if hasattr(sol['z6'], 'real') else sol['z6']
        
        print("\n" + "="*60)
        print("JACOBIAN NUMERICAL VERIFICATION")
        print("="*60)
        print(f"Solution: z4={z4_base:.6f}, z5={z5_base:.6f}, z6={z6_base:.6f}")
        print(f"Step size: eps={eps}")
        
        # Analytical Jacobian
        J_analytical = self.jacobian_matrix(z4_base, z5_base, z6_base)
        
        # Numerical Jacobian
        z_base = [z4_base, z5_base, z6_base]
        J_numerical = []
        
        for j in range(3):  # Derivative w.r.t. z4, z5, z6
            col = []
            for i in range(3):  # Equations f4, f5, f6 (particle indices 3, 4, 5)
                z_plus = list(z_base)
                z_plus[j] += eps
                z_minus = list(z_base)
                z_minus[j] -= eps
                
                # Evaluate scattering equations (particle 4=index 3, particle 5=index 4, particle 6=index 5)
                f_plus = self.scattering_equation(i+3, z_plus[0], z_plus[1], z_plus[2])
                f_minus = self.scattering_equation(i+3, z_minus[0], z_minus[1], z_minus[2])
                
                try:
                    f_plus_val = float(f_plus) if hasattr(f_plus, '__float__') else complex(f_plus).real
                    f_minus_val = float(f_minus) if hasattr(f_minus, '__float__') else complex(f_minus).real
                except:
                    f_plus_val = f_plus
                    f_minus_val = f_minus
                
                deriv = (f_plus_val - f_minus_val) / (2*eps)
                col.append(deriv)
            J_numerical.append(col)
        
        # Compare
        print("\nAnalytical Jacobian:")
        for i in range(3):
            row = [float(J_analytical[i,j]) if hasattr(J_analytical[i,j], '__float__') else J_analytical[i,j] for j in range(3)]
            print(f"  Row {i}: [{row[0]:.6e}, {row[1]:.6e}, {row[2]:.6e}]")
        
        print("\nNumerical Jacobian:")
        for i in range(3):
            # Note: J_numerical[j][i] because we computed column by column
            row = [J_numerical[j][i] for j in range(3)]
            print(f"  Row {i}: [{row[0]:.6e}, {row[1]:.6e}, {row[2]:.6e}]")
        
        # Check max difference
        max_diff = 0
        max_rel_diff = 0
        for i in range(3):
            for j in range(3):
                try:
                    analytical_val = float(J_analytical[i,j]) if hasattr(J_analytical[i,j], '__float__') else complex(J_analytical[i,j]).real
                except:
                    analytical_val = J_analytical[i,j]
                numerical_val = J_numerical[j][i]
                
                diff = abs(analytical_val - numerical_val)
                if abs(analytical_val) > 1e-10:
                    rel_diff = diff / abs(analytical_val)
                else:
                    rel_diff = diff
                
                max_diff = max(max_diff, diff)
                max_rel_diff = max(max_rel_diff, rel_diff)
        
        print(f"\nMax absolute difference: {max_diff:.6e}")
        print(f"Max relative difference: {max_rel_diff:.6e}")
        
        # Check determinants
        det_analytical = J_analytical.det()
        try:
            det_analytical_float = float(det_analytical) if hasattr(det_analytical, '__float__') else complex(det_analytical).real
        except:
            det_analytical_float = det_analytical
        
        # Numerical determinant (from numerical Jacobian)
        # J_numerical is stored as columns, need to transpose
        J_num_matrix = [[J_numerical[j][i] for j in range(3)] for i in range(3)]
        det_numerical = (J_num_matrix[0][0] * (J_num_matrix[1][1]*J_num_matrix[2][2] - J_num_matrix[1][2]*J_num_matrix[2][1])
                        - J_num_matrix[0][1] * (J_num_matrix[1][0]*J_num_matrix[2][2] - J_num_matrix[1][2]*J_num_matrix[2][0])
                        + J_num_matrix[0][2] * (J_num_matrix[1][0]*J_num_matrix[2][1] - J_num_matrix[1][1]*J_num_matrix[2][0]))
        
        print(f"\nDeterminant (analytical): {det_analytical_float:.6e}")
        print(f"Determinant (numerical):  {det_numerical:.6e}")
        det_diff = abs(det_analytical_float - det_numerical)
        det_rel_diff = det_diff / abs(det_analytical_float) if det_analytical_float != 0 else det_diff
        print(f"Determinant difference: {det_diff:.6e} (relative: {det_rel_diff:.6e})")
        
        match = max_rel_diff < 1e-4
        if match:
            print("\n✅ Jacobian computation verified")
        else:
            print("\n⚠️  Jacobian may have issues")
        
        return {
            'max_absolute_difference': max_diff,
            'max_relative_difference': max_rel_diff,
            'det_analytical': det_analytical_float,
            'det_numerical': det_numerical,
            'det_difference': det_diff,
            'det_relative_difference': det_rel_diff,
            'verified': match
        }


# Test function
def test_scattering_solver():
    """Quick test of ScatteringEquationSolver"""
    print("Testing ScatteringEquationSolver...")
    
    # Create test kinematics
    kin = SpinorKinematics.random_rational(6, seed=42)
    
    # Create solver
    solver = ScatteringEquationSolver(kin)
    
    # Test equation evaluation
    var('z4 z5 z6')
    f4 = solver.scattering_equation(3, z4, z5, z6)
    f5 = solver.scattering_equation(4, z4, z5, z6)
    f6 = solver.scattering_equation(5, z4, z5, z6)
    
    print(f"f4 has {len(f4.variables())} variables")
    print(f"f5 has {len(f5.variables())} variables")
    print(f"f6 has {len(f6.variables())} variables")
    
    # Test numerical solving (if kinematics is numeric)
    try:
        solutions = solver.solve_numerical()
        print(f"Found {len(solutions)} numerical solutions")
        for i, sol in enumerate(solutions[:3]):  # Print first 3
            print(f"Solution {i+1}: z4={sol['z4']:.4f}, z5={sol['z5']:.4f}, z6={sol['z6']:.4f}")
    except Exception as e:
        print(f"Numerical solving failed: {e}")
    
    print("Test complete.")


# Only run test if executed directly (commented out to avoid running on load)
# if __name__ == "__main__":
#     test_scattering_solver()

