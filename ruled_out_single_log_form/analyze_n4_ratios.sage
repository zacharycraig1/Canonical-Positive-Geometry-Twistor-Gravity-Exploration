#!/usr/bin/env sage
from sage.all import *

load('tests/test_n4_spinor_helicity.sage')

print("Analyzing n=4 ratios...")

ratios = []
for seed in range(20):
    try:
        lambdas, tildelambdas = sample_spinor_helicity_4pt(seed=seed)
        
        # Check angle brackets
        all_ok = True
        for i in range(4):
            for j in range(i+1, 4):
                if angle_bracket_spinor(lambdas[i], lambdas[j]) == 0:
                    all_ok = False
                    break
            if not all_ok:
                break
        if not all_ok:
            continue
        
        M4_klt = klt_4pt_spinor(lambdas, tildelambdas)
        M4_hodges = hodges_4pt_spinor(lambdas, tildelambdas)
        
        if M4_klt is None or M4_hodges is None or M4_hodges == 0:
            continue
        
        ratio = M4_klt / M4_hodges
        ratios.append((seed, ratio, abs(ratio)))
    except:
        continue

print(f"Collected {len(ratios)} ratios")

# Check if absolute values are constant
abs_ratios = [abs(r[1]) for r in ratios]
unique_abs = list(set(abs_ratios))
print(f"Unique absolute ratios: {len(unique_abs)}")

# Check if ratios are related by sign
signs = [1 if r[1] > 0 else -1 for r in ratios]
print(f"Sign pattern: {signs[:10]}")

# Check if ratios are constant up to a sign
if len(unique_abs) == 1:
    print(f"SUCCESS: Ratios are constant up to sign: |ratio| = {unique_abs[0]}")
elif len(unique_abs) <= 3:
    print(f"Few absolute values: {unique_abs}")
    # Try to see if they're related by simple factors
    for i, abs_r in enumerate(unique_abs[:3]):
        for j, abs_r2 in enumerate(unique_abs[:3]):
            if i != j and abs_r != 0 and abs_r2 != 0:
                factor = abs_r2 / abs_r
                print(f"  |r_{j}| / |r_{i}| = {factor}")

# Check individual components
print("\nFirst 3 samples:")
for i, (seed, ratio, abs_r) in enumerate(ratios[:3]):
    lambdas, tildelambdas = sample_spinor_helicity_4pt(seed=seed)
    M4_klt = klt_4pt_spinor(lambdas, tildelambdas)
    M4_hodges = hodges_4pt_spinor(lambdas, tildelambdas)
    print(f"\nSeed {seed}:")
    print(f"  KLT = {M4_klt}")
    print(f"  Hodges = {M4_hodges}")
    print(f"  Ratio = {ratio}")
    print(f"  |Ratio| = {abs_r}")

