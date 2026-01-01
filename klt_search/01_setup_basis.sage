#!/usr/bin/env sage
from sage.all import *
import itertools
import os

# =============================================================================
# STEP 1: DEFINE BASIS AND DIMENSIONS
# =============================================================================

def setup_basis():
    print("Step 1: Setting up KLT Ansatz Basis")
    
    # 1. Permutation Basis (S3)
    # Fixed legs: {1, 5, 6} (indices 0, 4, 5)
    # Permuted legs: {2, 3, 4} (indices 1, 2, 3)
    perm_indices = [1, 2, 3]
    perms = sorted(list(itertools.permutations(perm_indices)))
    print(f"Permutation Basis Size: {len(perms)} (alpha) x {len(perms)} (beta) = {len(perms)**2}")

    # 2. Monomial Basis (Degree 3 in 9 variables)
    # We use 9 independent Mandelstam invariants for n=6.
    # Basis Choice: s12, s23, s34, s45, s56, s61, s123, s234, s345
    # Note: This set might have linear dependencies? 
    # n(n-3)/2 = 9.
    # The set of planar adjacent 2-particle s_{i,i+1} is 6 variables.
    # The set of planar adjacent 3-particle s_{i,i+1,i+2} is 6 variables (s123=s456).
    # s123, s234, s345 are independent.
    # So {s12, s23, s34, s45, s56, s61, s123, s234, s345} is a candidate basis.
    # We will assume these 9 are our variables.
    
    num_vars = 9
    degree = 3
    
    # Generate exponent tuples for degree <= 3 (or exactly 3? KLT is homogenous degree 3 in momenta? s is momenta^2. 
    # KLT kernel S is order s^3? 
    # For n=4, S ~ s (dim 2). 
    # For n=5, S ~ s^2 (dim 4).
    # For n=6, S ~ s^3 (dim 6).
    # Yes, degree 3. It should be homogeneous of degree 3 if we consider dimensions.
    # But let's allow <= 3 to be safe, or just degree 3 if we are confident.
    # Let's allow all monomials up to degree 3.
    
    monomials = []
    # Use combinations_with_replacement to generate exponents summing to d
    from itertools import combinations_with_replacement
    
    # Map index 0..8 to the variable list
    basis_names = ["s12", "s23", "s34", "s45", "s56", "s61", "s123", "s234", "s345"]
    
    for d in range(degree + 1):
        for combo in combinations_with_replacement(range(num_vars), d):
            # Convert combo to exponent vector
            exponents = [0] * num_vars
            for idx in combo:
                exponents[idx] += 1
            monomials.append(tuple(exponents))
            
    print(f"Monomial Basis Size (Degree <= {degree}): {len(monomials)}")
    
    total_unknowns = (len(perms) ** 2) * len(monomials)
    print(f"Total Unknowns to Solve: {total_unknowns}")
    
    # Save to dictionary
    basis_data = {
        "perms": perms,
        "monomials": monomials,
        "basis_names": basis_names,
        "degree": degree,
        "num_vars": num_vars,
        "total_unknowns": total_unknowns
    }
    
    save(basis_data, "klt_search/basis_metadata.sobj")
    print("Saved basis metadata to klt_search/basis_metadata.sobj")

if __name__ == "__main__":
    setup_basis()









