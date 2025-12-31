from itertools import combinations

def get_unique_3pt_channels(n=6):
    """Returns list of unique 3-particle channels (pairs of complement sets)."""
    indices = set(range(n))
    seen_channels = set()
    channels_to_test = []
    
    for comb in combinations(indices, 3):
        c = tuple(sorted(comb))
        complement = tuple(sorted(list(indices - set(c))))
        
        if c in seen_channels or complement in seen_channels:
            continue
            
        seen_channels.add(c)
        seen_channels.add(complement)
        channels_to_test.append(c)
        
    return channels_to_test

def get_recommended_shift(channel, n=6):
    """
    Returns (shift_a, shift_b) such that shift_a in channel, shift_b out.
    """
    in_set = set(channel)
    out_set = set(range(n)) - in_set
    
    return list(in_set)[0], list(out_set)[0]



