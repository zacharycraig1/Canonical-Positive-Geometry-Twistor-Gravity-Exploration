
from sage.all import *
from itertools import permutations

def analyze_matrix_product():
    print("Analyzing M * S product structure...")
    
    data = load('klt_search/matrices_MS.sobj')
    M = data['M']
    S = data['S']
    basis = data['basis']
    
    P = M * S
    print("\nProduct M * S:")
    print(P)
    
    # Is it a scaled permutation matrix?
    # Or block diagonal?
    # Or just messy?
    
    # Normalize by first element if non-zero
    ref = P[0,0]
    if ref != 0:
        P_norm = P / ref
        print("\nNormalized P (approx):")
        # Print approximate values for easier reading
        for row in P_norm:
            print([float(x) for x in row])
            
    # Check if S is related to M^-1
    try:
        Minv = M.inverse()
        print("\nM_inv * S:")
        Prod2 = Minv * S
        # Normalize
        if Prod2[0,0] != 0:
            print("\nNormalized M_inv * S:")
            for row in (Prod2 / Prod2[0,0]):
                print([float(x) for x in row])
    except:
        print("M singular")

    # Try brute force permutation finding
    # Find P_pi such that P_pi * M * P_pi.T * S = c * I ? 
    # Or M * P_pi * S = c * I?
    # KLT kernel definition might use a permuted basis relative to Parke-Taylor.
    
    # Let's try to permute rows/cols of S to make M*S diagonal
    # S_new = P * S * P.T ? 
    # Or just permute indices of S.
    
    # M is indexed by [alpha | beta].
    # S is indexed by [alpha | beta].
    # If standard KLT formula assumes a different ordering of basis elements than our lexicographic one.
    
    # Try all 6! permutations of the basis applied to S
    # Actually, just permuting columns of P?
    # If M * S = D * Perm, then S = M^-1 * D * Perm.
    
    # Let's check if rows of P are proportional to unit vectors.
    # From output above: (1, 1, 1, 0, 0, 0)
    # This looks like mixing.
    
    # Check rank of S
    print(f"\nRank of S: {S.rank()}")
    if S.rank() < 6:
        print("S is singular!")
        
    # Inspect first few rows of S to see if columns are identical
    print("\nS matrix (first 2 rows):")
    print(S[:2])
    
    # Check if cols 0,1,2 are identical
    col0 = S.column(0)
    col1 = S.column(1)
    col2 = S.column(2)
    if col0 == col1: print("Col 0 == Col 1")
    if col0 == col2: print("Col 0 == Col 2")
    
    # If S is singular, my KLT implementation is broken or using redundant basis.
    
    pass

if __name__ == "__main__":
    analyze_matrix_product()

