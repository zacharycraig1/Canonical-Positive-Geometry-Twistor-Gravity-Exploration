# src/kinematic_associahedron/cell_decomposition.sage
"""
Cell Decomposition of Gravity Amplitude
========================================

The key hypothesis: The Hodges amplitude can be written as a sum over
"cells" (triangulations of the associahedron), similar to how:

    m(α,α) = Σ_triangulations 1/(X_1 × X_2 × X_3)

For gravity, we expect something like:

    M_6 = <12>^8 × Σ_triangulations N_T / (X_1 × X_2 × X_3)

where N_T is a numerator that depends on the triangulation.

If the gravity amplitude is a canonical form, then this cell structure
should emerge naturally from the geometry.

Strategy:
1. Compute Hodges at multiple kinematic points
2. For each triangulation T, compute 1/(X_1 × X_2 × X_3)
3. Try to express Hodges as linear combination of these terms
4. Find the coefficients (numerators)
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


def square_bracket(tilde_lambdas, i, j):
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def compute_planar_mandelstam(i, j, kin):
    """
    Compute the planar Mandelstam X_ij = s_{i,i+1,...,j-1}.
    
    For ordering (0,1,2,3,4,5), X_ij is the sum of s_ab for particles
    between position i and position j (exclusive of j).
    """
    n = 6
    
    # Indices in the channel
    if i < j:
        indices = list(range(i, j))
    else:
        indices = list(range(i, n)) + list(range(0, j))
    
    # Sum s_ab for all pairs in this set
    s_total = QQ(0)
    for idx_a in range(len(indices)):
        for idx_b in range(idx_a + 1, len(indices)):
            a = indices[idx_a]
            b = indices[idx_b]
            s_total += kin.s(a, b)
    
    return s_total


def get_triangulations():
    """
    Get all 14 triangulations of a hexagon.
    
    Each triangulation is a set of 3 non-crossing diagonals.
    A diagonal (i,j) with j > i+1 and (i,j) ≠ (0,5).
    """
    n = 6
    
    # All diagonals (non-edges)
    diagonals = []
    for i in range(n):
        for j in range(i+2, n):
            if not (i == 0 and j == n-1):
                diagonals.append((i, j))
    
    # Check if two diagonals cross
    def crosses(d1, d2):
        a, b = d1
        c, d = d2
        if a > b: a, b = b, a
        if c > d: c, d = d, c
        return (a < c < b < d) or (c < a < d < b)
    
    # Find all triangulations (maximal non-crossing sets)
    triangulations = []
    
    def backtrack(current, start_idx):
        if len(current) == n - 3:  # 3 diagonals for hexagon
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
    """
    Compute 1 / (X_1 × X_2 × X_3) for a triangulation.
    
    Each diagonal (i,j) corresponds to X_ij = s_{i,...,j-1}.
    """
    denom = QQ(1)
    
    for (i, j) in triangulation:
        X_ij = compute_planar_mandelstam(i, j, kin)
        if X_ij == 0:
            return None  # On a pole
        denom *= X_ij
    
    return QQ(1) / denom


def analyze_cell_structure(num_samples=5):
    """
    Analyze if Hodges can be written as sum over cells.
    """
    print("\n" + "="*70)
    print("CELL DECOMPOSITION ANALYSIS")
    print("="*70)
    
    triangulations = get_triangulations()
    print(f"\nNumber of triangulations: {len(triangulations)}")
    
    print("\nTriangulations:")
    for i, T in enumerate(triangulations):
        print(f"  T_{i}: {T}")
    
    # For each kinematic point, compute:
    # 1. Hodges amplitude
    # 2. Each triangulation term 1/(X_1 X_2 X_3)
    # 3. See if Hodges = Σ c_T × triangulation_term
    
    all_data = []
    
    for seed in range(42, 42 + num_samples):
        kin = SpinorKinematics.random_rational(n=6, seed=seed)
        
        M_6, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
        
        if M_6 is None:
            print(f"\nSeed {seed}: Hodges failed [{status}]")
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
            print(f"\nSeed {seed}: On a pole")
            continue
        
        # Extract helicity factor
        ang_12 = angle_bracket(kin.lambdas, 0, 1)
        h_factor = ang_12**8
        
        # Reduced amplitude (remove helicity)
        if h_factor == 0:
            continue
        M_reduced = M_6 / h_factor
        
        all_data.append({
            'seed': seed,
            'M_6': M_6,
            'M_reduced': M_reduced,
            'h_factor': h_factor,
            'terms': terms
        })
        
        print(f"\nSeed {seed}:")
        print(f"  M_6 = {float(M_6):.6e}")
        print(f"  <12>^8 = {float(h_factor):.6e}")
        print(f"  M_6/<12>^8 = {float(M_reduced):.6e}")
        print(f"  Triangulation terms (1/X_1X_2X_3):")
        for i, term in enumerate(terms):
            print(f"    T_{i}: {float(term):.6e}")
    
    return all_data, triangulations


def try_linear_decomposition(all_data, triangulations):
    """
    Try to find coefficients c_T such that:
    M_reduced = Σ c_T × triangulation_term
    
    If the coefficients are CONSTANT across all kinematic points,
    then we've found the cell decomposition.
    """
    print("\n" + "="*70)
    print("TRYING LINEAR DECOMPOSITION")
    print("="*70)
    
    if len(all_data) < 2:
        print("Need at least 2 data points")
        return
    
    n_triangulations = len(triangulations)
    
    # Build a system: for each data point,
    # M_reduced = c_0 × term_0 + c_1 × term_1 + ... + c_13 × term_13
    
    # If we have more data points than unknowns (14), we can check consistency
    
    print(f"\nNumber of data points: {len(all_data)}")
    print(f"Number of unknowns (coefficients): {n_triangulations}")
    
    # Build matrix A where A[i][j] = term_j for data point i
    # and vector b where b[i] = M_reduced for data point i
    
    A = matrix(QQ, len(all_data), n_triangulations)
    b = vector(QQ, len(all_data))
    
    for i, data in enumerate(all_data):
        for j, term in enumerate(data['terms']):
            A[i, j] = term
        b[i] = data['M_reduced']
    
    print(f"\nMatrix A shape: {A.nrows()} × {A.ncols()}")
    print(f"Rank of A: {A.rank()}")
    
    if A.rank() < n_triangulations:
        print("Matrix is rank-deficient, cannot solve uniquely")
    
    # Try to solve Ax = b
    try:
        # If overdetermined, use least squares
        if A.nrows() > A.ncols():
            # Check if b is in the column space of A
            Ab = A.augment(b.column())
            if Ab.rank() > A.rank():
                print("System is inconsistent (no solution)")
                
                # Find the least-squares residual
                AtA = A.transpose() * A
                Atb = A.transpose() * b
                
                if AtA.rank() == n_triangulations:
                    c = AtA.solve_right(Atb)
                    residual = A * c - b
                    print(f"Least-squares solution found")
                    print(f"Residual norm: {float(residual.norm()):.6e}")
                else:
                    print("Cannot compute least-squares solution")
                    return
            else:
                # Consistent system, find any solution
                c = A.solve_right(b)
                print("Exact solution found!")
        else:
            c = A.solve_right(b)
            print("Solution found!")
        
        print(f"\nCoefficients c_T:")
        for i in range(n_triangulations):
            print(f"  c_{i} = {c[i]} = {float(c[i]):.6e}")
        
        # Verify the solution
        print(f"\nVerification:")
        for i, data in enumerate(all_data):
            predicted = sum(c[j] * data['terms'][j] for j in range(n_triangulations))
            actual = data['M_reduced']
            error = abs(float(predicted - actual))
            print(f"  Seed {data['seed']}: predicted={float(predicted):.6e}, actual={float(actual):.6e}, error={error:.6e}")
            
    except Exception as e:
        print(f"Could not solve: {e}")


def check_bcj_structure():
    """
    Check if the amplitude has BCJ (color-kinematics) structure.
    
    In BCJ, the amplitude can be written as:
    M = Σ_diagrams n_i^2 / D_i
    
    where n_i are numerators satisfying Jacobi identities,
    and D_i are products of propagators.
    """
    print("\n" + "="*70)
    print("CHECKING BCJ STRUCTURE")
    print("="*70)
    
    print("\nBCJ (Bern-Carrasco-Johansson) duality states:")
    print("  Gravity = (YM numerator)² / propagators")
    print("")
    print("For MHV gravity:")
    print("  M_n = Σ_diagrams n_i² / (s_1 × s_2 × ... × s_{n-3})")
    print("")
    print("The numerators n_i satisfy Jacobi identities:")
    print("  n_s + n_t + n_u = 0")
    print("")
    print("This suggests the gravity amplitude is a SUM over cubic diagrams,")
    print("which are exactly the triangulations of the associahedron!")
    print("")
    print("If BCJ holds, then:")
    print("  M_6 = <12>^8 × Σ_{14 triangulations} n_T² / (X_1 × X_2 × X_3)")
    print("")
    print("where n_T are YM numerators for each triangulation.")


def main():
    """Run the full cell decomposition analysis."""
    # Analyze cell structure
    all_data, triangulations = analyze_cell_structure(num_samples=10)
    
    # Try linear decomposition
    if len(all_data) >= 2:
        try_linear_decomposition(all_data, triangulations)
    
    # Explain BCJ connection
    check_bcj_structure()
    
    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print("""
If the linear decomposition gives CONSTANT coefficients, then:
  M_6 = <12>^8 × Σ_T c_T / (X_1 × X_2 × X_3)

These coefficients c_T would be the "cell weights" of the positive geometry.

If the coefficients VARY with kinematics, then:
- The gravity amplitude is NOT a simple sum over associahedron cells
- The geometry is more complex (e.g., fiber bundle, or different polytope)
- BCJ numerators may depend on kinematics in a specific way
""")


if __name__ == "__main__":
    main()

