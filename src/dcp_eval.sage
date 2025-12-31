from sage.all import *
from itertools import combinations

try:
    from src.kinematics_map import spinors_to_channels, spinors_to_sij
except ImportError:
    # If run via load(), these might already be available or need loading
    # Ideally load() in the caller handles dependencies.
    pass

# =============================================================================
# DCP CHANNEL GENERATION (1-based to match 54.sage)
# =============================================================================

def canonical_three_subset(S, n=6):
    S_set = set(S)
    comp = set(range(1, n+1)) - S_set
    # Sort to compare
    t1 = tuple(sorted(S_set))
    t2 = tuple(sorted(comp))
    return min(t1, t2)

def generate_channels_n6():
    """
    Generate sorted list of channels C for n=6 (1-based indices).
    Matches 54.sage logic.
    """
    C = []
    # 2-particle channels
    for ij in combinations(range(1,7), 2):
        C.append(tuple(sorted(ij)))
    # 3-particle channels
    for S in combinations(range(1,7), 3):
        C.append(canonical_three_subset(S, 6))
        
    return sorted(list(set(C)))

def generate_triples(C):
    """
    Generate all triples of channel indices (i,j,k) with i<j<k.
    This is the basis for the naive 3-form space.
    """
    m = len(C)
    triples = []
    for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                triples.append((i, j, k))
    return triples

# =============================================================================
# EVALUATION
# =============================================================================

def eval_form_on_channels(candidate_vec, channels_dict, C, triples):
    """
    Evaluate form on channel values.
    
    Args:
        candidate_vec: List/Vector of coefficients c_{ijk} for each triple.
        channels_dict: Dict mapping channel tuple (1-based) to value (QQ).
        C: List of channel tuples (indexed 0..m-1).
        triples: List of (i,j,k) indices into C.
        
    Returns:
        Value sum c_{ijk} / (s_i * s_j * s_k)
    """
    result = QQ(0)
    
    for idx, (i, j, k) in enumerate(triples):
        if idx >= len(candidate_vec):
            break
        
        coeff = candidate_vec[idx]
        if coeff == 0:
            continue
            
        # Get channel tuples
        ch_i = C[i]
        ch_j = C[j]
        ch_k = C[k]
        
        # Get values
        # channels_dict keys must match C entries
        si = channels_dict.get(ch_i)
        sj = channels_dict.get(ch_j)
        sk = channels_dict.get(ch_k)
        
        if si is None or sj is None or sk is None:
            # Missing channel value? Should not happen if dict is complete
            continue
            
        if si == 0 or sj == 0 or sk == 0:
            # Singular point - skip or return infinity?
            # For evaluation on generic points, this shouldn't happen.
            continue
            
        term = coeff / (si * sj * sk)
        result += term
        
    return result

def eval_form_on_spinors(candidate_vec, lambdas, tilde_lambdas, C=None, triples=None):
    """
    Evaluate form on spinors.
    1. Compute s_ij and s_S (0-based).
    2. Convert to 1-based map.
    3. Call eval_form_on_channels.
    """
    if C is None:
        C = generate_channels_n6()
    if triples is None:
        triples = generate_triples(C)
        
    # 1. Compute 0-based channels
    # spinors_to_channels returns dict with 0-based keys
    channels_0 = spinors_to_channels(lambdas, tilde_lambdas)
    
    # 2. Convert to 1-based dict matching C format
    channels_1 = {}
    n = 6
    
    # Helper to canonicalize 1-based key
    def to_canon_1(key_0):
        # Convert 0-based tuple to 1-based set
        s_1 = set([x + 1 for x in key_0])
        
        # Only canonicalize against complement for n/2 size (here 3)
        if len(s_1) == 3:
            comp = set(range(1, n+1)) - s_1
            t1 = tuple(sorted(s_1))
            t2 = tuple(sorted(comp))
            return min(t1, t2)
        else:
            # For 2-particle channels, always use the set itself (sorted)
            return tuple(sorted(s_1))

    for k0, val in channels_0.items():
        k1 = to_canon_1(k0)
        # Store. Note: multiple 0-based keys might map to same 1-based canonical key 
        # (e.g. complementary sets). Value should be same.
        channels_1[k1] = val
        
    # 3. Evaluate
    return eval_form_on_channels(candidate_vec, channels_1, C, triples)

