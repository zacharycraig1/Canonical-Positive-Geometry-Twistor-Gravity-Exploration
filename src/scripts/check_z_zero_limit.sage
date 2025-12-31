import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_polynomial

def check_z_zero_limit():
    print("Checking Limit z_3_4 -> 0 in Forest Polynomial F(z)...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Full Polynomial
    F = get_forest_polynomial(n, roots)
    
    # 2. Limit z_3_4 -> 0
    # This effectively removes any forest containing edge (3,4)
    # The result should be the forest polynomial of the graph K_6 \ { (3,4) }
    
    # Find variable name for (3,4)
    # edge vars are z_i_j with i < j
    var_name = "z_3_4"
    R = F.parent()
    z34 = R(var_name)
    
    F_limit = F.subs({z34: 0})
    
    print(f"Original F terms: {len(F.monomials())}")
    print(f"Limit F terms: {len(F_limit.monomials())}")
    
    # Verify: The removed terms are exactly those with z_3_4
    # i.e. Forests containing edge (3,4).
    # If we interpret (3,4) as a factorization channel?
    # If [3,4] -> 0, does the amplitude factorize?
    # Usually [3,4]->0 means "anti-collinear".
    # Factorization involves contracting the edge?
    
    # 3. Geometric Interpretation
    # The facet z_3_4 >= 0 corresponds to the polynomial F_limit.
    # This polynomial enumerates forests on the graph with edge (3,4) DELETED.
    # Is this a product of lower point objects?
    # If deleting (3,4) disconnects the graph? No, K_6 is highly connected.
    # But for trees, deleting an edge splits it into two components.
    # For FORESTS (rooted), deleting an edge might split a component?
    # But the condition is "each component has exactly one root".
    # If we delete (3,4) from K_6, we just lose some valid forests.
    # Does F_limit factorize?
    
    factors = F_limit.factor()
    print(f"Factors object: {factors}")
    
    # Check if non-trivial factorization (more than 1 factor, ignoring unit)
    # factor() returns Factorization object which acts like list of (poly, exp)
    if len(factors) > 1 or (len(factors) == 1 and factors[0][1] > 1):
         print("F_limit FACTORIZES! This suggests physical factorization.")
    else:
         print("F_limit does not factorize. So z_3_4 -> 0 is not a factorization channel.")
        
    
    # 4. Check Infinity Limit (Coefficient of z_3_4)
    # Physically z -> infinity corresponds to <ij> -> 0 (holomorphic collinear)
    print("\nChecking Limit z_3_4 -> infinity (Coefficient of z_3_4)...")
    F_infty = F.coefficient({z34: 1})
    print(f"Terms in coefficient: {len(F_infty.monomials())}")
    
    factors_inf = F_infty.factor()
    print(f"Factors: {factors_inf}")
    
    if len(factors_inf) > 1 or (len(factors_inf) == 1 and factors_inf[0][1] > 1):
         print("F_infty FACTORIZES! This confirms physical factorization at <34> -> 0.")
    else:
         print("F_infty does not factorize.")
         
    print("\nConclusion: The facet z_ij = 0 corresponds to deleting an edge.")
    
if __name__ == "__main__":
    check_z_zero_limit()

