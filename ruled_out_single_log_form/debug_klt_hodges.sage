#!/usr/bin/env sage
# =============================================================================
# DEBUG: Why is Hodges returning 0?
# =============================================================================

from sage.all import *
import numpy as np
from itertools import combinations, permutations

# Load the KLT proof functions
load('correct_klt_proof.sage')

# Test on a single moment-curve point
Z = sample_positive_Z_moment_curve(n=6, seed=0)
twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)

print("="*70)
print("DEBUGGING SINGLE POINT")
print("="*70)
print(f"Domain OK: {twistor.domain_ok}")
print(f"Domain reason: {twistor.domain_reason}")

# Check some brackets
print("\nSome angle brackets:")
for i in range(6):
    j = (i + 1) % 6
    ang = twistor.get_angle(i, j)
    print(f"  <{i},{j}> = {ang}")

print("\nSome four-brackets:")
for ijkl in [(0,1,2,3), (1,2,3,4), (2,3,4,5)]:
    bracket = twistor.get_four_bracket(*ijkl)
    print(f"  <{ijkl}> = {bracket}")

# Compute Hodges
print("\nComputing Hodges...")
H_result = hodges_6pt_mhv(twistor)
H = H_result[0] if isinstance(H_result, tuple) else H_result
H_reason = H_result[1] if isinstance(H_result, tuple) else "ok"

print(f"Hodges = {H}")
print(f"H_reason = {H_reason}")

# Check the Phi matrix
print("\nComputing Phi matrix manually...")
n = 6
indices = [1, 2, 3, 4]
d = len(indices)

Phi = matrix(QQ, d, d)

for ii, i in enumerate(indices):
    for jj, j in enumerate(indices):
        if ii == jj:
            diag_sum = QQ(0)
            for k in range(n):
                if k in [i, 0, 5]:
                    continue
                ik_sq = twistor.get_square(i, k)
                ik_ang = twistor.get_angle(i, k)
                i0_ang = twistor.get_angle(i, 0)
                i5_ang = twistor.get_angle(i, 5)
                k0_ang = twistor.get_angle(k, 0)
                k5_ang = twistor.get_angle(k, 5)
                if ik_ang == 0 or i0_ang == 0 or i5_ang == 0:
                    continue
                contrib = ik_sq * k0_ang * k5_ang / (ik_ang * i0_ang * i5_ang)
                diag_sum -= contrib
            Phi[ii, jj] = diag_sum
        else:
            ij_ang = twistor.get_angle(i, j)
            ij_sq = twistor.get_square(i, j)
            Phi[ii, jj] = ij_sq / ij_ang

print(f"Phi matrix:\n{Phi}")
print(f"det(Phi) = {Phi.det()}")

denom = QQ(1)
for i in range(n):
    j = (i + 1) % n
    bracket = twistor.get_angle(i, j)
    denom *= bracket

print(f"Denominator (product of <i i+1>) = {denom}")
print(f"Hodges = det(Phi) / denom = {Phi.det() / denom}")

# Compute KLT
print("\nComputing KLT...")
A_result = gravity_6pt_mhv_klt(twistor)
A = A_result[0] if isinstance(A_result, tuple) else A_result
A_reason = A_result[1] if isinstance(A_result, tuple) else "ok"

print(f"KLT = {A}")
print(f"A_reason = {A_reason}")

if H != 0 and A is not None:
    ratio = A / H
    print(f"\nRatio A/H = {ratio}")
    print(f"Difference = {A - H}")

