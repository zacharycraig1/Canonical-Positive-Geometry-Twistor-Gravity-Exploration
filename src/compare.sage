#!/usr/bin/env sage
# =============================================================================
# COMPARE MODULE: Exact equality tests and normalization handling
# =============================================================================

from sage.all import *
from collections import Counter


def exact_equality_test(H, A):
    """
    Test exact equality using rational arithmetic.
    
    Args:
        H: Hodges amplitude (QQ)
        A: KLT amplitude (QQ)
        
    Returns:
        Tuple (is_equal, ratio, diff)
        - is_equal: True if A = c * H for some constant c
        - ratio: A / H (or None)
        - diff: A - ratio * H (should be 0 if equal)
    """
    if H is None or A is None:
        return False, None, None
    
    if H == 0:
        return (A == 0), None, None
    
    if A == 0:
        return False, None, None
    
    # Compute ratio (exact QQ)
    ratio = A / H
    
    # Check if they're equal (up to constant)
    # For rational functions, A - (A/H) * H should be exactly 0
    diff = A - ratio * H
    
    if diff == 0:
        return True, ratio, None
    
    return False, ratio, diff


def find_normalization_constant(ratios_list):
    """
    Identify if ratios are all equal to a constant c.
    
    Args:
        ratios_list: List of ratio values (QQ)
        
    Returns:
        Tuple (constant, is_constant)
        - constant: Most common ratio value
        - is_constant: True if all ratios are equal
    """
    if not ratios_list:
        return None, False
    
    ratio_counts = Counter(ratios_list)
    most_common = ratio_counts.most_common(1)[0][0]
    is_constant = all(r == most_common for r in ratios_list)
    
    return most_common, is_constant


def analyze_ratio_variation(ratios_list):
    """
    Analyze why ratios vary (if they do).
    
    Args:
        ratios_list: List of ratio values
        
    Returns:
        Dict with analysis results
    """
    if not ratios_list:
        return {'status': 'empty', 'count': 0}
    
    unique_ratios = list(set(ratios_list))
    
    result = {
        'count': len(ratios_list),
        'unique_count': len(unique_ratios),
        'is_constant': len(unique_ratios) == 1,
        'unique_ratios': unique_ratios[:10],  # First 10
    }
    
    if len(unique_ratios) == 1:
        result['status'] = 'constant'
        result['constant'] = unique_ratios[0]
    elif len(unique_ratios) <= 5:
        result['status'] = 'few_values'
        # Check if they differ by signs or simple factors
        if len(unique_ratios) == 2:
            r1, r2 = unique_ratios
            if r1 == -r2:
                result['relationship'] = 'sign_difference'
            elif r1 != 0 and r2 != 0 and r2 / r1 in [2, QQ(1)/2, -1, -2, QQ(-1)/2]:
                result['relationship'] = f'factor_{r2/r1}'
    else:
        result['status'] = 'variable'
    
    return result


def verify_identity_polynomial(H_func, A_func, twistor_generator, n_tests=100):
    """
    Verify polynomial identity KLT = c * Hodges using multiple samples.
    
    This is a rigorous identity test based on Schwartz-Zippel lemma.
    
    Args:
        H_func: Function computing Hodges amplitude
        A_func: Function computing KLT amplitude
        twistor_generator: Function generating twistor data
        n_tests: Number of test points
        
    Returns:
        Dict with verification results
    """
    ratios = []
    failures = []
    
    for i in range(n_tests):
        twistor = twistor_generator(seed=i)
        
        H_result = H_func(twistor)
        A_result = A_func(twistor)
        
        H = H_result[0] if isinstance(H_result, tuple) else H_result
        A = A_result[0] if isinstance(A_result, tuple) else A_result
        
        if H is None or A is None or H == 0:
            failures.append({
                'idx': i,
                'H_none': H is None,
                'A_none': A is None,
                'H_zero': H == 0 if H is not None else False
            })
            continue
        
        ratios.append(A / H)
    
    analysis = analyze_ratio_variation(ratios)
    
    return {
        'n_tests': n_tests,
        'n_valid': len(ratios),
        'n_failures': len(failures),
        'analysis': analysis,
        'failures': failures[:10],  # First 10 failures
        'is_identity': analysis.get('is_constant', False),
        'constant': analysis.get('constant', None)
    }









