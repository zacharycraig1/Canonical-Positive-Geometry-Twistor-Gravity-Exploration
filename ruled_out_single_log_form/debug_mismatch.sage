#!/usr/bin/env sage
from sage.all import *

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')

print("Finding mismatch in random integer sampling...")

# Find the failing random integer case
for seed in range(50):
    twistor = MomentumTwistor(n=6, seed=seed, check_domain=True)
    if not twistor.domain_ok:
        print(f"Seed {seed}: domain violation")
        continue
    H_result = hodges_6pt_mhv(twistor)
    A_result = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
    H = H_result[0] if isinstance(H_result, tuple) else H_result
    A = A_result[0] if isinstance(A_result, tuple) else A_result
    if H is None or A is None:
        print(f"Seed {seed}: H={H}, A={A}")
        continue
    is_equal, ratio, diff = exact_equality_test(H, A)
    if not is_equal:
        print(f"MISMATCH at seed {seed}")
        print(f"H = {H}")
        print(f"A = {A}")
        print(f"ratio = {ratio}")
        print(f"diff = {diff}")
        break
    else:
        print(f"Seed {seed}: OK, ratio={ratio}")

