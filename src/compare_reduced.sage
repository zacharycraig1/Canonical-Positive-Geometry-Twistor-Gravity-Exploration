#!/usr/bin/env sage
# =============================================================================
# REDUCED-ONLY COMPARE MODULE
# =============================================================================
# Dedicated comparison function for proof mode:
#   R(Z) := M_6^{KLT}(Z) / \bar M_6^{Hodges}(Z)
# 
# NO wrappers, NO cyclic prefactors, NO det' conversions.
# Pure reduced-to-reduced comparison.
# =============================================================================

from sage.all import *


def compare_reduced_only(klt_result, hodges_reduced_result):
    """
    Compare KLT gravity to Hodges reduced amplitude.
    
    Computes: R(Z) = M_6^{KLT}(Z) / \bar M_6^{Hodges}(Z)
    
    Args:
        klt_result: Tuple (amplitude, reason) from gravity_6pt_mhv_klt
        hodges_reduced_result: Tuple (amplitude, reason) from hodges_6pt_mhv_reduced
        
    Returns:
        Tuple (ratio, is_valid, reason)
        - ratio: R(Z) = KLT / HodgesReduced (or None if invalid)
        - is_valid: True if both amplitudes computed successfully
        - reason: "ok" or error description
    """
    # Extract values
    klt_val = klt_result[0] if isinstance(klt_result, tuple) else klt_result
    hodges_val = hodges_reduced_result[0] if isinstance(hodges_reduced_result, tuple) else hodges_reduced_result
    
    # Check for None
    if klt_val is None:
        return (None, False, "klt_none")
    if hodges_val is None:
        return (None, False, "hodges_none")
    
    # Check for zero
    if hodges_val == 0:
        return (None, False, "hodges_zero")
    
    # Compute ratio
    ratio = klt_val / hodges_val
    
    return (ratio, True, "ok")


def test_reduced_ratio_constant(twistor_generator, n_samples=50):
    """
    Test if reduced ratio is constant across samples.
    
    Args:
        twistor_generator: Function(seed) -> MomentumTwistor
        n_samples: Number of samples to test
        
    Returns:
        Tuple (is_constant, constant_value, unique_count, ratio_list)
    """
    ratios = []
    
    for seed in range(n_samples):
        twistor = twistor_generator(seed)
        
        if not twistor.domain_ok:
            continue
        
        # Import here to avoid circular imports
        from src.hodges import hodges_6pt_mhv_reduced
        from src.klt import gravity_6pt_mhv_klt, mandelstam_invariant
        
        H_red = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        
        ratio, is_valid, reason = compare_reduced_only(A_klt, H_red)
        
        if is_valid:
            ratios.append(ratio)
    
    if not ratios:
        return (False, None, 0, [])
    
    unique_ratios = list(set(ratios))
    is_constant = len(unique_ratios) == 1
    
    return (is_constant, unique_ratios[0] if is_constant else None, len(unique_ratios), ratios)










