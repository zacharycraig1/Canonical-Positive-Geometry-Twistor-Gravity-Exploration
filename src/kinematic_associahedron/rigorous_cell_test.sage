# src/kinematic_associahedron/rigorous_cell_test.sage
"""
RIGOROUS TEST: Cell Decomposition Hypothesis

Hypothesis: The 6-point MHV gravity amplitude can be written as:
    M_6 / <12>^8 = Σ_{i=0}^{13} c_i / (X_{a_i} × X_{b_i} × X_{c_i})

where c_i are CONSTANT coefficients (independent of kinematics).

To test this properly, we need:
1. MORE data points than unknowns (at least 15, preferably 25+)
2. Check that the matrix A has full rank (rank = 14)
3. Check that the residual is essentially zero (not just "small")
4. CRITICAL: Check that coefficients are CONSISTENT across subsets of data

This is a proper statistical test, not just curve fitting.
"""

from sage.all import *
from itertools import permutations, combinations
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


def angle_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def compute_planar_mandelstam(i, j, kin):
    """Compute the planar Mandelstam X_ij = s_{i,i+1,...,j-1}."""
    n = 6
    if i < j:
        indices = list(range(i, j))
    else:
        indices = list(range(i, n)) + list(range(0, j))
    
    s_total = QQ(0)
    for idx_a in range(len(indices)):
        for idx_b in range(idx_a + 1, len(indices)):
            a = indices[idx_a]
            b = indices[idx_b]
            s_total += kin.s(a, b)
    
    return s_total


def get_triangulations():
    """Get all 14 triangulations of a hexagon."""
    n = 6
    diagonals = []
    for i in range(n):
        for j in range(i+2, n):
            if not (i == 0 and j == n-1):
                diagonals.append((i, j))
    
    def crosses(d1, d2):
        a, b = d1
        c, d = d2
        if a > b: a, b = b, a
        if c > d: c, d = d, c
        return (a < c < b < d) or (c < a < d < b)
    
    triangulations = []
    
    def backtrack(current, start_idx):
        if len(current) == n - 3:
            triangulations.append(tuple(sorted(current)))
            return
        
        for i in range(start_idx, len(diagonals)):
            d = diagonals[i]
            is_valid = all(not crosses(d, existing) for existing in current)
            if is_valid:
                backtrack(current + [d], i + 1)
    
    backtrack([], 0)
    return triangulations


def compute_triangulation_term(triangulation, kin):
    """Compute 1 / (X_1 × X_2 × X_3) for a triangulation."""
    denom = QQ(1)
    for (i, j) in triangulation:
        X_ij = compute_planar_mandelstam(i, j, kin)
        if X_ij == 0:
            return None
        denom *= X_ij
    return QQ(1) / denom


def collect_data(num_samples, seed_start=42):
    """Collect kinematic data for regression analysis."""
    triangulations = get_triangulations()
    n_triangulations = len(triangulations)
    
    print(f"\nCollecting {num_samples} kinematic samples...")
    print(f"Number of triangulations (unknowns): {n_triangulations}")
    
    all_data = []
    
    for seed in range(seed_start, seed_start + num_samples * 2):  # Try more seeds to get enough valid ones
        if len(all_data) >= num_samples:
            break
            
        try:
            kin = SpinorKinematics.random_rational(n=6, seed=seed)
            
            M_6, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
            
            if M_6 is None:
                continue
            
            # Compute triangulation terms
            terms = []
            valid = True
            for T in triangulations:
                term = compute_triangulation_term(T, kin)
                if term is None:
                    valid = False
                    break
                terms.append(term)
            
            if not valid:
                continue
            
            # Extract helicity factor
            ang_12 = angle_bracket(kin.lambdas, 0, 1)
            h_factor = ang_12**8
            
            if h_factor == 0:
                continue
                
            M_reduced = M_6 / h_factor
            
            all_data.append({
                'seed': seed,
                'M_reduced': M_reduced,
                'terms': terms
            })
            
        except Exception as e:
            continue
    
    print(f"Collected {len(all_data)} valid samples")
    return all_data, triangulations


def rigorous_regression_test(all_data, triangulations):
    """
    Perform rigorous regression test.
    
    Key checks:
    1. Is the system overdetermined? (n_samples > n_unknowns)
    2. Is the design matrix full rank?
    3. Is the residual truly zero?
    4. Are coefficients stable across data subsets?
    """
    print("\n" + "="*70)
    print("RIGOROUS REGRESSION TEST")
    print("="*70)
    
    n_samples = len(all_data)
    n_unknowns = len(triangulations)
    
    print(f"\nSamples: {n_samples}")
    print(f"Unknowns: {n_unknowns}")
    
    if n_samples < n_unknowns:
        print(f"\n❌ UNDERDETERMINED SYSTEM: {n_samples} < {n_unknowns}")
        print("Cannot uniquely determine coefficients!")
        print("Any 'solution' found is just ONE particular solution, not THE solution.")
        return None, False
    
    # Build design matrix A and target vector b
    A = matrix(QQ, n_samples, n_unknowns)
    b = vector(QQ, n_samples)
    
    for i, data in enumerate(all_data):
        for j, term in enumerate(data['terms']):
            A[i, j] = term
        b[i] = data['M_reduced']
    
    rank_A = A.rank()
    print(f"\nMatrix A shape: {A.nrows()} × {A.ncols()}")
    print(f"Rank of A: {rank_A}")
    
    if rank_A < n_unknowns:
        print(f"\n❌ RANK DEFICIENT: rank {rank_A} < {n_unknowns}")
        print("Columns are linearly dependent - coefficients not uniquely determined!")
        return None, False
    
    print(f"\n✓ System is overdetermined and full rank")
    
    # Solve via least squares: minimize ||Ax - b||^2
    AtA = A.transpose() * A
    Atb = A.transpose() * b
    
    try:
        c = AtA.solve_right(Atb)
    except Exception as e:
        print(f"\n❌ Could not solve: {e}")
        return None, False
    
    # Compute residual
    residual = A * c - b
    residual_norm = sqrt(sum(r**2 for r in residual))
    max_abs_residual = max(abs(float(r)) for r in residual)
    
    print(f"\nResidual L2 norm: {float(residual_norm):.6e}")
    print(f"Max absolute residual: {max_abs_residual:.6e}")
    
    # Check if residual is essentially zero (machine precision for rationals)
    is_exact = (residual_norm == 0)
    
    if is_exact:
        print("\n✓ EXACT FIT: Residual is exactly zero!")
    else:
        print(f"\n❌ NON-ZERO RESIDUAL: The hypothesis is FALSIFIED!")
        print("The amplitude cannot be written as a sum over triangulations with constant coefficients.")
        return c, False
    
    # Print coefficients
    print("\nCoefficients:")
    for i in range(n_unknowns):
        c_val = c[i]
        if c_val == 0:
            print(f"  c_{i} = 0 (exactly)")
        else:
            print(f"  c_{i} = {float(c_val):.6e}")
    
    # Count zero and non-zero coefficients
    n_zero = sum(1 for cv in c if cv == 0)
    n_nonzero = n_unknowns - n_zero
    print(f"\nNon-zero coefficients: {n_nonzero}")
    print(f"Zero coefficients: {n_zero}")
    
    return c, True


def cross_validation_test(all_data, triangulations):
    """
    Cross-validation: solve on subsets and check consistency.
    
    If the decomposition is true, coefficients from ANY overdetermined 
    subset should match coefficients from any other subset.
    """
    print("\n" + "="*70)
    print("CROSS-VALIDATION TEST")
    print("="*70)
    
    n_samples = len(all_data)
    n_unknowns = len(triangulations)
    
    if n_samples < n_unknowns + 5:
        print(f"Need at least {n_unknowns + 5} samples for cross-validation")
        return False
    
    # Split data into two halves
    half = n_samples // 2
    
    if half < n_unknowns:
        print(f"Each half has only {half} samples, need {n_unknowns}")
        return False
    
    data_1 = all_data[:half]
    data_2 = all_data[half:]
    
    print(f"\nSplit: {len(data_1)} + {len(data_2)} samples")
    
    def solve_subset(data):
        A = matrix(QQ, len(data), n_unknowns)
        b = vector(QQ, len(data))
        for i, d in enumerate(data):
            for j, term in enumerate(d['terms']):
                A[i, j] = term
            b[i] = d['M_reduced']
        
        if A.rank() < n_unknowns:
            return None
        
        AtA = A.transpose() * A
        Atb = A.transpose() * b
        return AtA.solve_right(Atb)
    
    c1 = solve_subset(data_1)
    c2 = solve_subset(data_2)
    
    if c1 is None or c2 is None:
        print("One or both subsets are rank-deficient")
        return False
    
    # Check if coefficients match
    diff = c1 - c2
    max_diff = max(abs(float(d)) for d in diff)
    
    print(f"\nMax coefficient difference between subsets: {max_diff:.6e}")
    
    if max_diff == 0:
        print("✓ EXACT MATCH: Coefficients are identical across subsets!")
        return True
    elif max_diff < 1e-10:
        print("✓ Coefficients match to high precision")
        return True
    else:
        print("❌ COEFFICIENTS DIFFER: Decomposition hypothesis is FALSIFIED!")
        print("\nCoefficient comparison:")
        for i in range(n_unknowns):
            print(f"  c_{i}: subset1={float(c1[i]):.6e}, subset2={float(c2[i]):.6e}, diff={float(diff[i]):.6e}")
        return False


def verify_on_new_data(coefficients, triangulations, num_tests=10, seed_start=200):
    """
    Final verification: use computed coefficients to predict NEW data points.
    """
    print("\n" + "="*70)
    print("PREDICTION TEST ON NEW DATA")
    print("="*70)
    
    if coefficients is None:
        print("No coefficients to test")
        return False
    
    successes = 0
    attempts = 0
    
    for seed in range(seed_start, seed_start + num_tests * 2):
        if attempts >= num_tests:
            break
            
        try:
            kin = SpinorKinematics.random_rational(n=6, seed=seed)
            M_6, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
            
            if M_6 is None:
                continue
            
            terms = []
            valid = True
            for T in triangulations:
                term = compute_triangulation_term(T, kin)
                if term is None:
                    valid = False
                    break
                terms.append(term)
            
            if not valid:
                continue
            
            ang_12 = angle_bracket(kin.lambdas, 0, 1)
            h_factor = ang_12**8
            if h_factor == 0:
                continue
            
            M_reduced_actual = M_6 / h_factor
            M_reduced_predicted = sum(coefficients[j] * terms[j] for j in range(len(triangulations)))
            
            diff = M_reduced_actual - M_reduced_predicted
            
            attempts += 1
            
            if diff == 0:
                print(f"  Seed {seed}: ✓ EXACT MATCH")
                successes += 1
            else:
                rel_error = abs(float(diff / M_reduced_actual)) if M_reduced_actual != 0 else float('inf')
                print(f"  Seed {seed}: ❌ MISMATCH (rel error: {rel_error:.6e})")
                
        except Exception as e:
            continue
    
    print(f"\nPrediction accuracy: {successes}/{attempts}")
    
    if successes == attempts and attempts > 0:
        print("✓ ALL PREDICTIONS EXACT!")
        return True
    else:
        print("❌ Some predictions failed")
        return False


def main():
    print("="*70)
    print("RIGOROUS CELL DECOMPOSITION TEST")
    print("="*70)
    print("""
This test determines whether the 6-point MHV gravity amplitude
can be expressed as a sum over the 14 triangulations of the 
kinematic associahedron with CONSTANT coefficients:

    M_6 / <12>^8 = Σ_{i=0}^{13} c_i / (X_{a_i} × X_{b_i} × X_{c_i})

If true, this would identify the positive geometry for gravity.
""")
    
    # Collect 25 samples (well over the 14 unknowns)
    all_data, triangulations = collect_data(num_samples=30, seed_start=42)
    
    if len(all_data) < 14:
        print(f"\n❌ Only {len(all_data)} valid samples - cannot run test")
        return
    
    # Test 1: Rigorous regression
    coefficients, is_exact = rigorous_regression_test(all_data, triangulations)
    
    if not is_exact:
        print("\n" + "="*70)
        print("CONCLUSION: HYPOTHESIS FALSIFIED")
        print("="*70)
        print("""
The amplitude CANNOT be written as a sum over triangulations 
with constant coefficients. The residual is non-zero.
""")
        return
    
    # Test 2: Cross-validation
    cv_passed = cross_validation_test(all_data, triangulations)
    
    if not cv_passed:
        print("\n" + "="*70)
        print("CONCLUSION: HYPOTHESIS FALSIFIED")
        print("="*70)
        print("""
Coefficients are NOT consistent across data subsets.
The decomposition does not hold.
""")
        return
    
    # Test 3: Prediction on new data
    pred_passed = verify_on_new_data(coefficients, triangulations, num_tests=10, seed_start=300)
    
    print("\n" + "="*70)
    if is_exact and cv_passed and pred_passed:
        print("CONCLUSION: HYPOTHESIS VERIFIED!")
        print("="*70)
        print("""
The 6-point MHV gravity amplitude CAN be expressed as a sum
over triangulations of the kinematic associahedron!

    M_6 = <12>^8 × Σ_{T} c_T / (X_1 × X_2 × X_3)

This identifies the positive geometry for gravity as the
kinematic associahedron with BCJ-like weights.
""")
    else:
        print("CONCLUSION: HYPOTHESIS NOT FULLY VERIFIED")
        print("="*70)
        print("Some tests passed, some failed. More investigation needed.")


if __name__ == "__main__":
    main()

