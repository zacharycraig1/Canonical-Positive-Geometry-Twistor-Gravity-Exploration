#!/usr/bin/env sage
# =============================================================================
# DEBUG: BCFW CELL COMPUTATION
# =============================================================================
# Test BCFW formula on a known good point to see what's happening
# =============================================================================

from sage.all import *
import numpy as np
from itertools import combinations

# Load the working Hodges implementation
load('correct_amplituhedron_hodges.sage')

# Create a test point that we know works
twistor = MomentumTwistor(n=6, seed=42)

print("Testing on point with seed=42")
print(f"Z[0] = {twistor.Z[0]}")
print(f"Z[1] = {twistor.Z[1]}")

# Test Hodges
H = hodges_6pt_mhv(twistor)
print(f"\nHodges = {H}")

# Test BCFW
A, cells = amplituhedron_from_bcfw_cells(twistor)
print(f"BCFW sum = {A}")
print(f"Number of cells = {len(cells)}")

if A is not None and H is not None:
    if H != 0:
        ratio = A / H
        print(f"\nRatio A/H = {ratio}")
        print(f"Difference = {A - H}")
    else:
        print("\nH is zero")
else:
    print(f"\nA is None: {A is None}")
    print(f"H is None: {H is None}")

# Check some brackets
print("\nChecking some brackets:")
print(f"<0,1> = {twistor.get_angle(0, 1)}")
print(f"<1,2> = {twistor.get_angle(1, 2)}")
print(f"<0,1,2,3> = {twistor.get_four_bracket(0, 1, 2, 3)}")

