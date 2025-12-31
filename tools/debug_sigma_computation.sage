#!/usr/bin/env sage
# =============================================================================
# DEBUG SIGMA COMPUTATION
# =============================================================================
# Manually verify sigma computation for a few cases
# =============================================================================

from sage.all import *
load('src/hodges_sigma.sage')

def manual_sigma_check(I, J, n):
    """Manually compute sigma to verify correctness."""
    print(f"\nComputing sigma for I={I}, J={J}, n={n}")
    
    # Convert to 1-indexed
    I_1 = [i + 1 for i in I]
    J_1 = [j + 1 for j in J]
    
    all_labels = list(range(1, n + 1))
    restI = sorted([x for x in all_labels if x not in I_1])
    restJ = sorted([x for x in all_labels if x not in J_1])
    
    A = I_1 + restI
    B = J_1 + restJ
    
    print(f"  I (1-indexed): {I_1}")
    print(f"  J (1-indexed): {J_1}")
    print(f"  restI: {restI}")
    print(f"  restJ: {restJ}")
    print(f"  A = {A}")
    print(f"  B = {B}")
    
    # Build position map
    posB = {val: idx for idx, val in enumerate(B)}
    p_list = [posB[val] for val in A]
    
    print(f"  Permutation (positions in B): {p_list}")
    
    # Compute sign
    sigma_func = hodges_sigma(I, J, n)
    
    # Manual sign computation
    inversions = 0
    for i in range(len(p_list)):
        for j in range(i + 1, len(p_list)):
            if p_list[i] > p_list[j]:
                inversions += 1
    
    sigma_manual = 1 if inversions % 2 == 0 else -1
    
    print(f"  Sigma (function): {sigma_func}")
    print(f"  Sigma (manual): {sigma_manual}")
    print(f"  Match: {sigma_func == sigma_manual}")
    
    return sigma_func

# Test cases
print("="*70)
print("DEBUGGING SIGMA COMPUTATION")
print("="*70)

manual_sigma_check([0,1,2], [3,4,5], 6)
manual_sigma_check([0,1,2], [0,1,2], 6)
manual_sigma_check([0,1,3], [2,4,5], 6)






