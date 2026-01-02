# Debug: Why is the residual so large?
# Let's check if the basic math is correct

from sage.all import *
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical

def angle_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]

def compute_planar_mandelstam(i, j, kin):
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

# Get a single kinematic point
kin = SpinorKinematics.random_rational(n=6, seed=42)
M_6, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))

print("=== Debug: Single Point Analysis ===\n")

print(f"M_6 = {float(M_6):.6e}")

ang_12 = angle_bracket(kin.lambdas, 0, 1)
h_factor = ang_12**8
M_reduced = M_6 / h_factor

print(f"<12>^8 = {float(h_factor):.6e}")
print(f"M_reduced = M_6 / <12>^8 = {float(M_reduced):.6e}")

triangulations = get_triangulations()
print(f"\nTriangulations: {len(triangulations)}")

print("\nTriangulation terms:")
for i, T in enumerate(triangulations):
    Xs = []
    product = QQ(1)
    for (a, b) in T:
        X = compute_planar_mandelstam(a, b, kin)
        Xs.append(X)
        product *= X
    term = QQ(1) / product
    print(f"  T_{i}: {T}")
    print(f"       X's: {[float(x) for x in Xs]}")
    print(f"       1/(X1*X2*X3) = {float(term):.6e}")

# Check order of magnitude
print(f"\n\nScale comparison:")
print(f"  M_reduced ~ {float(M_reduced):.2e}")
print(f"  Triangulation terms ~ {float(term):.2e}")
print(f"  Ratio ~ {float(M_reduced / term):.2e}")

# For the gravity amplitude, the BCJ claim is that there should be 
# kinematic-dependent NUMERATORS n_T, not constant coefficients
print("\n\nConclusion:")
print("If M_reduced / term is NOT constant across kinematics,")
print("then the simple constant-coefficient hypothesis is FALSE.")
print("The BCJ picture says: M = Σ n_T^2 / D_T")
print("where n_T are kinematic-dependent BCJ numerators.")

