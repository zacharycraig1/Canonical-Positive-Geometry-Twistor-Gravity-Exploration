#!/usr/bin/env sage
from sage.all import *

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')

print("Verifying constant ratio with 200 samples...")

ratios = []
for seed in range(200):
    Z = sample_positive_Z_moment_curve(n=6, seed=seed)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    if not twistor.domain_ok:
        continue
    H = hodges_6pt_mhv_reduced(twistor)
    A = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
    H_val = H[0] if isinstance(H, tuple) else H
    A_val = A[0] if isinstance(A, tuple) else A
    if H_val is None or A_val is None or H_val == 0:
        continue
    ratios.append(A_val / H_val)

unique = list(set(ratios))
print(f'Total samples: {len(ratios)}')
print(f'Unique ratios: {len(unique)}')
if len(unique) == 1:
    print(f'CONSTANT RATIO: {unique[0]}')
    print(f'Ratio (float): {float(unique[0]):.15e}')
    # Try to simplify
    try:
        simplified = unique[0].simplify_rational()
        print(f'Simplified: {simplified}')
    except:
        pass
else:
    print(f'First 3 unique: {unique[:3]}')

