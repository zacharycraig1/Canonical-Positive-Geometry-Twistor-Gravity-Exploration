#!/usr/bin/env sage
# =============================================================================
# Positivity Options Testing for 6-Point MHV Gravity
# =============================================================================
# Systematically tests the four positivity conditions from the directive:
# B1: All entries Psi_ij > 0
# B2: Pf'(Psi) > 0
# B3: All principal minors of Psi_reduced > 0 (totally positive)
# B4: The matrix (i*Psi) is positive semi-definite

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/psi_matrix.sage")

class PositivityOptionTester:
    """
    Tests all four positivity options (B1-B4) to determine which defines R6.
    """
    
    def __init__(self, kinematics):
        """
        Initialize with kinematics.
        
        Args:
            kinematics: SpinorKinematics object (n=6)
        """
        self.kin = kinematics
        
    def test_option_b1(self, z4, z5, z6):
        """
        Option B1: All (non-gauge-fixed) entries Psi_ij > 0
        
        From directive lines 396-401:
        For indices {4,5,6}, requires:
        - [45]²/(z_4 - z_5) > 0  =>  z_4 > z_5
        - [46]²/(z_4 - z_6) > 0  =>  z_4 > z_6
        - [56]²/(z_5 - z_6) > 0  =>  z_5 > z_6
        
        Combined: z_4 > z_5 > z_6 > 0
        
        Args:
            z4, z5, z6: worldsheet coordinates (numeric)
            
        Returns:
            True if point satisfies Option B1
        """
        # Basic ordering check
        if not (z4 > z5 > z6 > 0):
            return False
        
        # For MHV with particles 1,2 negative, 3-6 positive:
        # The square brackets [ij]² are generically positive for generic kinematics
        # So the sign depends only on (z_i - z_j)
        
        # The above ordering ensures all entries among {4,5,6} are positive
        
        # Mixed entries (between 1,2 and 4,5,6) have sign depending on kinematics
        # For now, we only enforce ordering constraints
        
        return True
    
    def test_option_b2(self, z4, z5, z6):
        """
        Option B2: Pf'(Psi) > 0
        
        Evaluates the reduced Pfaffian at the given point and checks sign.
        
        Args:
            z4, z5, z6: worldsheet coordinates (numeric)
            
        Returns:
            (is_positive, pfaffian_value) tuple
        """
        # Create Psi matrix at this point
        psi = PsiMatrixMHV(self.kin, z4, z5, z6)
        
        # Compute reduced Pfaffian
        try:
            pf_prime = psi.reduced_pfaffian_standard((0, 1))
            if pf_prime == 0:
                # Try alternative deletion
                pf_prime = psi.reduced_pfaffian_delete_3_6()
        except:
            try:
                pf_prime = psi.reduced_pfaffian_delete_3_6()
            except:
                return (None, None)  # Computation failed
        
        # Convert to float for comparison
        try:
            pf_val = float(pf_prime) if hasattr(pf_prime, '__float__') else complex(pf_prime).real
        except:
            pf_val = pf_prime
        
        is_positive = (pf_val > 0)
        
        return (is_positive, pf_val)
    
    def test_option_b3(self, z4, z5, z6):
        """
        Option B3: All principal minors of Psi_reduced > 0 (totally positive)
        
        For a 4x4 reduced matrix (after deleting rows/cols 1,2), check:
        - All 1x1 minors (diagonal entries) > 0
        - All 2x2 principal minors > 0
        - All 3x3 principal minors > 0
        - The 4x4 determinant > 0
        
        Args:
            z4, z5, z6: worldsheet coordinates (numeric)
            
        Returns:
            (is_totally_positive, details) tuple
        """
        # Create Psi matrix at this point
        psi = PsiMatrixMHV(self.kin, z4, z5, z6)
        
        # Get reduced matrix (delete rows/cols 1,2 -> indices {3,4,5,6} remain)
        # But with z3=infinity, row/col 3 is all zeros
        # So use deletion (1,6) or (2,5) to get non-degenerate matrix
        
        try:
            # Delete particles 1,6 (indices 0,5)
            M_red = psi.reduced_matrix([0, 5])
        except:
            return (None, "Matrix construction failed")
        
        if M_red.nrows() != 4:
            return (None, f"Expected 4x4 matrix, got {M_red.nrows()}x{M_red.ncols()}")
        
        # Check all principal minors
        details = {
            '1x1_minors': [],
            '2x2_minors': [],
            '3x3_minors': [],
            '4x4_det': None,
            'all_positive': True
        }
        
        # 1x1 principal minors (diagonal entries)
        # For antisymmetric matrix, diagonal is 0, so this is trivial
        # Instead check off-diagonal entries in upper triangle
        for i in range(4):
            for j in range(i+1, 4):
                entry = M_red[i, j]
                try:
                    entry_val = float(entry) if hasattr(entry, '__float__') else complex(entry).real
                except:
                    entry_val = entry
                details['1x1_minors'].append(entry_val)
        
        # 2x2 principal minors
        from itertools import combinations
        for idx_pair in combinations(range(4), 2):
            i, j = idx_pair
            minor_2x2 = M_red[[i,j], [i,j]].det()
            try:
                minor_val = float(minor_2x2) if hasattr(minor_2x2, '__float__') else complex(minor_2x2).real
            except:
                minor_val = minor_2x2
            details['2x2_minors'].append(minor_val)
            if minor_val <= 0:
                details['all_positive'] = False
        
        # 3x3 principal minors
        for idx_triple in combinations(range(4), 3):
            i, j, k = idx_triple
            minor_3x3 = M_red[[i,j,k], [i,j,k]].det()
            try:
                minor_val = float(minor_3x3) if hasattr(minor_3x3, '__float__') else complex(minor_3x3).real
            except:
                minor_val = minor_3x3
            details['3x3_minors'].append(minor_val)
            if minor_val <= 0:
                details['all_positive'] = False
        
        # 4x4 determinant
        det_4x4 = M_red.det()
        try:
            det_val = float(det_4x4) if hasattr(det_4x4, '__float__') else complex(det_4x4).real
        except:
            det_val = det_4x4
        details['4x4_det'] = det_val
        if det_val <= 0:
            details['all_positive'] = False
        
        return (details['all_positive'], details)
    
    def test_option_b4(self, z4, z5, z6):
        """
        Option B4: The matrix (i*Psi) is positive semi-definite
        
        For antisymmetric matrix Psi, i*Psi is Hermitian.
        Check if all eigenvalues >= 0.
        
        Args:
            z4, z5, z6: worldsheet coordinates (numeric)
            
        Returns:
            (is_psd, eigenvalues) tuple
        """
        # Create Psi matrix at this point
        psi = PsiMatrixMHV(self.kin, z4, z5, z6)
        
        # Get reduced matrix
        try:
            M_red = psi.reduced_matrix([0, 5])  # Delete particles 1,6
        except:
            return (None, None)
        
        # For antisymmetric real matrix A, i*A is Hermitian
        # Eigenvalues of i*A are purely imaginary multiples of A's eigenvalues
        # For PSD, we need all eigenvalues of i*A to be real and >= 0
        
        # Actually, for real antisymmetric A:
        # - Eigenvalues come in pairs ±iλ where λ is real
        # - i*A has eigenvalues that are real: ±λ
        # For i*A to be PSD, we'd need all eigenvalues real and non-negative
        # But antisymmetric matrices have imaginary eigenvalues, so i*A has real eigenvalues
        
        try:
            # Compute eigenvalues of M_red (antisymmetric)
            # Then eigenvalues of i*M_red
            eigenvals_Psi = M_red.eigenvalues()
            
            # For antisymmetric matrix, eigenvalues are purely imaginary or zero
            # Eigenvalues of i*Psi are real (i times imaginary = real)
            eigenvals_iPsi = [complex(i * ev).real for ev in eigenvals_Psi]
            
            # Check if all are non-negative
            is_psd = all(ev >= -1e-10 for ev in eigenvals_iPsi)  # Small tolerance for numerical error
            
            return (is_psd, eigenvals_iPsi)
        except:
            return (None, None)
    
    def test_all_options(self, z4, z5, z6, verbose=False):
        """
        Test all four options at a given point.
        
        Args:
            z4, z5, z6: worldsheet coordinates
            verbose: Print detailed results
            
        Returns:
            dict with results for each option
        """
        results = {}
        
        # Test B1
        results['B1'] = self.test_option_b1(z4, z5, z6)
        
        # Test B2
        b2_positive, b2_value = self.test_option_b2(z4, z5, z6)
        results['B2'] = {
            'is_positive': b2_positive,
            'pfaffian_value': b2_value
        }
        
        # Test B3
        b3_tp, b3_details = self.test_option_b3(z4, z5, z6)
        results['B3'] = {
            'is_totally_positive': b3_tp,
            'details': b3_details
        }
        
        # Test B4
        b4_psd, b4_eigenvals = self.test_option_b4(z4, z5, z6)
        results['B4'] = {
            'is_psd': b4_psd,
            'eigenvalues': b4_eigenvals
        }
        
        if verbose:
            print(f"\nTesting point (z4={z4:.4f}, z5={z5:.4f}, z6={z6:.4f}):")
            print(f"  B1 (Entry positivity):      {results['B1']}")
            pf_val = results['B2']['pfaffian_value']
            pf_str = f"{pf_val:.6e}" if pf_val is not None else "None"
            print(f"  B2 (Pfaffian positivity):   {results['B2']['is_positive']} (Pf'={pf_str})")
            print(f"  B3 (Total positivity):      {results['B3']['is_totally_positive']}")
            print(f"  B4 (PSD i*Psi):            {results['B4']['is_psd']}")
        
        return results
    
    def scan_region(self, option_name, z_range=(0.1, 3.0), num_points=10):
        """
        Scan a region of (z4, z5, z6) space to understand option's structure.
        
        Args:
            option_name: 'B1', 'B2', 'B3', or 'B4'
            z_range: (min, max) for z values
            num_points: number of points per dimension
            
        Returns:
            dict with statistics
        """
        z_min, z_max = z_range
        z_vals = [z_min + (z_max - z_min) * i / (num_points - 1) for i in range(num_points)]
        
        count_satisfies = 0
        count_total = 0
        
        for z4 in z_vals:
            for z5 in z_vals:
                for z6 in z_vals:
                    # Basic ordering for sampling
                    if not (z4 > z5 > z6 > 0):
                        continue
                    
                    count_total += 1
                    
                    if option_name == 'B1':
                        satisfies = self.test_option_b1(z4, z5, z6)
                    elif option_name == 'B2':
                        satisfies, _ = self.test_option_b2(z4, z5, z6)
                    elif option_name == 'B3':
                        satisfies, _ = self.test_option_b3(z4, z5, z6)
                    elif option_name == 'B4':
                        satisfies, _ = self.test_option_b4(z4, z5, z6)
                    else:
                        raise ValueError(f"Unknown option: {option_name}")
                    
                    if satisfies:
                        count_satisfies += 1
        
        fraction = count_satisfies / count_total if count_total > 0 else 0
        
        return {
            'option': option_name,
            'total_points_tested': count_total,
            'points_satisfying': count_satisfies,
            'fraction': fraction
        }


# Test function
def test_positivity_options():
    """Quick test of PositivityOptionTester"""
    print("Testing PositivityOptionTester...")
    
    from src.kinematics.spinors import SpinorKinematics
    
    # Create test kinematics
    kin = SpinorKinematics.random_rational(6, seed=42)
    
    # Create tester
    tester = PositivityOptionTester(kin)
    
    # Test at a few points
    test_points = [
        (2.0, 1.5, 1.0),
        (1.5, 1.0, 0.5),
        (0.8, 0.5, 0.2),
    ]
    
    for z4, z5, z6 in test_points:
        results = tester.test_all_options(z4, z5, z6, verbose=True)
    
    # Scan regions
    print("\n" + "="*60)
    print("Scanning regions for each option:")
    for option in ['B1', 'B2', 'B3', 'B4']:
        stats = tester.scan_region(option, num_points=5)
        frac = float(stats['fraction']) * 100
        print(f"\n{option}: {stats['points_satisfying']}/{stats['total_points_tested']} points satisfy ({frac:.1f}%)")
    
    print("\nTest complete.")


# Only run test if executed directly
if __name__ == "__main__":
    test_positivity_options()

