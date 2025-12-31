#!/usr/bin/env sage
# =============================================================================
# QUICK TEST: Amplituhedron for 6-point MHV Gravity
# =============================================================================
# Simplified test to verify the amplituhedron approach works
# =============================================================================

from sage.all import *
import numpy as np

print("\n" + "="*70)
print("TESTING AMPLITUHEDRON APPROACH")
print("="*70)

# Generate simple momentum twistor kinematics
np.random.seed(42)
n = 6

# Create twistors
Z = []
for i in range(n):
    z = vector(QQ, [QQ(np.random.randint(-5, 6)), 
                     QQ(np.random.randint(-5, 6)),
                     QQ(np.random.randint(-5, 6)),
                     QQ(np.random.randint(-5, 6))])
    while all(x == 0 for x in z):
        z = vector(QQ, [QQ(np.random.randint(-5, 6)), 
                         QQ(np.random.randint(-5, 6)),
                         QQ(np.random.randint(-5, 6)),
                         QQ(np.random.randint(-5, 6))])
    Z.append(z)

print(f"\nGenerated {n} momentum twistors")

# Compute angle brackets
angle = {}
for i in range(n):
    for j in range(n):
        angle[(i,j)] = Z[i][0] * Z[j][1] - Z[i][1] * Z[j][0]

print("Computed angle brackets")

# Compute 4-brackets
four_bracket = {}
for ijkl in combinations(range(n), 4):
    i, j, k, l = sorted(ijkl)
    M = matrix(QQ, [Z[i], Z[j], Z[k], Z[l]])
    four_bracket[ijkl] = M.det()

print(f"Computed {len(four_bracket)} four-brackets")

# Test Hodges formula (simplified)
print("\nTesting Hodges formula...")
# Build Phi matrix for particles 2,3,4,5
indices = [1, 2, 3, 4]
Phi = matrix(QQ, 4, 4)

for ii, i in enumerate(indices):
    for jj, j in enumerate(indices):
        if ii == jj:
            Phi[ii, jj] = QQ(1)  # Simplified diagonal
        else:
            ij_angle = angle[(i, j)]
            if ij_angle == 0:
                print(f"  Singular: <{i}{j}> = 0")
                Phi[ii, jj] = QQ(0)
            else:
                Phi[ii, jj] = QQ(1) / ij_angle  # Simplified

try:
    det_Phi = Phi.det()
    print(f"  det(Phi) = {det_Phi}")
except:
    print("  Failed to compute det(Phi)")

# Test amplituhedron cells
print("\nTesting amplituhedron cells...")
cells = []
for i in range(n):
    for j in range(i+2, n):
        if j != (i + n - 1) % n:
            cells.append((i, j))

print(f"  Found {len(cells)} cells: {cells}")

# Compute cell forms
print("\nComputing cell forms...")
for cell in cells[:5]:  # Just first 5
    i, j = cell
    ip1 = (i + 1) % n
    jp1 = (j + 1) % n
    
    # Get 4-bracket
    indices = tuple(sorted([i, ip1, j, jp1]))
    if indices in four_bracket:
        form_val = four_bracket[indices]
        print(f"  Cell {cell}: <{i}{ip1}{j}{jp1}> = {form_val}")
    else:
        print(f"  Cell {cell}: bracket not found")

print("\n" + "="*70)
print("TEST COMPLETE")
print("="*70)

