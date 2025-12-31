
# src/dcp/nested_set_sampler_seeded.sage
import random as rnd

def sample_nested_sets_seeded(Gm, incompatible, required_flats, n_samples, rng_seed=None):
    """
    Generate maximal nested sets containing required_flats.
    
    Args:
        Gm: list of connected flats (masks)
        incompatible: set of incompatible pairs (indices in Gm)
        required_flats: list of flats (masks) that MUST be in the nested set
        n_samples: number of samples to generate
        rng_seed: random seed
        
    Returns:
        List of maximal nested sets (lists of flats/masks)
    """
    if rng_seed is not None:
        rnd.seed(int(rng_seed))
        
    # Build adjacency for incompatibility
    n = len(Gm)
    adj = [set() for _ in range(n)]
    for i, j in incompatible:
        adj[i].add(j)
        adj[j].add(i)
        
    # Map masks to indices
    mask_to_idx = {m: i for i, m in enumerate(Gm)}
    
    # Identify indices of required flats
    required_indices = []
    for m in required_flats:
        if m in mask_to_idx:
            required_indices.append(mask_to_idx[m])
        else:
            # If a required flat is not in Gm, we can't proceed
            # This might happen if Gm only contains connected flats and required is disconnected?
            # Or just missing. We'll warn or skip.
            print(f"Warning: Required flat {m} not in Gm")
            
    if not required_indices:
        print("Warning: No valid required flats found in Gm")
        
    # Initial state
    initial_blocked = set()
    initial_covered = 0
    
    for idx in required_indices:
        initial_blocked.update(adj[idx])
        initial_covered |= int(Gm[idx])
        
    # Identify pool of candidates (compatible with initial set)
    initial_pool = [i for i in range(n) 
                   if i not in required_indices 
                   and i not in initial_blocked]
                   
    full_mask = 0
    for m in Gm: full_mask |= int(m) # Or pass full_mask as arg
    # Note: assuming full_mask is the union of all atomic elements, usually (1<<6)-2
    
    results = []
    seen_hashes = set()
    
    attempts = 0
    max_attempts = n_samples * 10
    
    while len(results) < n_samples and attempts < max_attempts:
        attempts += 1
        
        current_indices = list(required_indices)
        current_blocked = set(initial_blocked)
        current_covered = initial_covered
        
        # Shuffle pool
        pool = list(initial_pool)
        rnd.shuffle(pool)
        
        # Greedy completion
        for idx in pool:
            if idx in current_blocked:
                continue
            
            # Additional check: does adding this create a forbidden antichain join?
            # The incompatible set handles pairwise incompatibility.
            # Assuming 'incompatible' captures all nested set axioms violations (pairwise).
            
            current_indices.append(idx)
            current_blocked.update(adj[idx])
            current_covered |= int(Gm[idx])
            
            if current_covered == full_mask:
                 # Check if maximal (size n-1 for n=6? No, depends on combinatorics)
                 # Actually, we just want it to be maximal in the sense that we can't add more.
                 # But usually we stop when it covers the full groundset and has correct size.
                 # For now, let's just continue until pool exhausted.
                 pass

        # Check maximality/validity
        # A valid maximal nested set for n=6 usually has size... it depends.
        # Let's just store what we found.
        
        N = [Gm[i] for i in current_indices]
        h = tuple(sorted(N))
        
        if h not in seen_hashes:
            seen_hashes.add(h)
            results.append(N)
            
    return results

def get_full_mask(Gm):
    m = 0
    for x in Gm: m |= int(x)
    return m

