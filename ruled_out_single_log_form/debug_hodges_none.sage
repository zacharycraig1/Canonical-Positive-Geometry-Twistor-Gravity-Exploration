#!/usr/bin/env sage
# Quick debug: Why is Hodges returning None?

load('proper_proof_amplituhedron_hodges.sage')

# Test on a single moment-curve point
Z = sample_positive_Z_moment_curve(n=6, seed=0)
twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)

print("="*70)
print("DEBUG: Why is Hodges returning None?")
print("="*70)
print(f"Domain OK: {twistor.domain_ok}")
print(f"Domain reason: {twistor.domain_reason}")

# Check some brackets
print("\nSome angle brackets:")
for i in range(6):
    j = (i + 1) % 6
    ang = twistor.get_angle(i, j)
    print(f"  <{i},{j}> = {ang}")

# Check normalization factor brackets
print("\nNormalization factor brackets:")
print(f"  <0,1> = {twistor.get_angle(0, 1)}")
print(f"  <1,2> = {twistor.get_angle(1, 2)}")
print(f"  <2,0> = {twistor.get_angle(2, 0)}")

# Compute Hodges
print("\nComputing Hodges...")
H_result = hodges_6pt_mhv(twistor)
H = H_result[0] if isinstance(H_result, tuple) else H_result
H_reason = H_result[1] if isinstance(H_result, tuple) else "ok"

print(f"Hodges = {H}")
print(f"H_reason = {H_reason}")

if H is None:
    print(f"\nERROR: Hodges returned None with reason: {H_reason}")

