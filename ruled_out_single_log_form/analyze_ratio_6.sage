#!/usr/bin/env sage
from sage.all import *

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')

print("Analyzing 6-ratio pattern...")

# Collect ratios with seed info
ratio_data = []
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
    ratio = A_val / H_val
    ratio_data.append({
        'seed': seed,
        'ratio': ratio,
        'ratio_str': str(ratio)
    })

# Group by ratio
from collections import defaultdict
ratio_groups = defaultdict(list)
for d in ratio_data:
    ratio_groups[d['ratio_str']].append(d['seed'])

print(f"\nFound {len(ratio_groups)} unique ratios:")
for i, (ratio_str, seeds) in enumerate(sorted(ratio_groups.items())[:6]):
    print(f"\nRatio {i+1}: {len(seeds)} occurrences")
    print(f"  Value: {ratio_str[:80]}...")
    print(f"  Seeds (first 10): {seeds[:10]}")
    print(f"  Seed pattern: mod 10 = {[s % 10 for s in seeds[:10]]}")
    print(f"  Seed pattern: mod 100 = {[s % 100 for s in seeds[:10]]}")

# Check if ratios are related
ratios_list = [QQ(r) for r in ratio_groups.keys()]
if len(ratios_list) >= 2:
    r0, r1 = ratios_list[0], ratios_list[1]
    print(f"\nChecking ratio relationships:")
    print(f"  r1/r0 = {r1/r0}")
    print(f"  r0/r1 = {r0/r1}")

