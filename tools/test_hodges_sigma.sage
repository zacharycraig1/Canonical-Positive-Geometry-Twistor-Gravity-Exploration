#!/usr/bin/env sage
# =============================================================================
# TEST HODGES SIGMA FUNCTION
# =============================================================================

from sage.all import *
load('src/hodges_sigma.sage')

print("Testing hodges_sigma function...")

# Test case 1: Standard deletion (0,1,2) and (3,4,5)
I1 = [0, 1, 2]
J1 = [3, 4, 5]
sigma1 = hodges_sigma(I1, J1, 6)
print(f"sigma([0,1,2], [3,4,5], 6) = {sigma1}")

# Test case 2: Symmetric deletion (0,1,2) and (0,1,2)
I2 = [0, 1, 2]
J2 = [0, 1, 2]
sigma2 = hodges_sigma(I2, J2, 6)
print(f"sigma([0,1,2], [0,1,2], 6) = {sigma2}")

# Test case 3: Mixed deletion (0,1,3) and (2,4,5)
I3 = [0, 1, 3]
J3 = [2, 4, 5]
sigma3 = hodges_sigma(I3, J3, 6)
print(f"sigma([0,1,3], [2,4,5], 6) = {sigma3}")

# Test case 4: Another mixed (0,2,4) and (1,3,5)
I4 = [0, 2, 4]
J4 = [1, 3, 5]
sigma4 = hodges_sigma(I4, J4, 6)
print(f"sigma([0,2,4], [1,3,5], 6) = {sigma4}")

print("\nExpected: All should be ±1")









