# dcp_flats.sage
# Assumes dcp_common.sage is loaded

def get_connected_flats(M6):
    cp = cache_path(CACHE_FLATS_FILE)
    if CACHE_PRECOMPUTE and os.path.exists(cp):
        log("Loading flats from cache...")
        data = load(cp)
        return data['Gm'], data['full_mask']
    
    log("Building connected flats...")
    t0 = time.time()
    full_set = set(M6.groundset())
    full_mask = set_to_mask(full_set)
    Gm, seen = [], set()
    
    for k in range(int(M6.rank()) + 1):
        for F in M6.flats(k):
            if not F: continue
            Fm = set_to_mask(F)
            if Fm in (PY0, full_mask) or Fm in seen: continue
            seen.add(Fm)
            if M6.delete(list(full_set - set(F))).is_connected():
                Gm.append(Fm)
    
    Gm.sort(key=lambda x: (popcount(x), int(x)))
    log(f"Built {len(Gm)} flats in {time.time()-t0:.1f}s")
    if CACHE_PRECOMPUTE:
        save({'Gm': Gm, 'full_mask': int(full_mask)}, cp)
    return Gm, int(full_mask)

def get_incompatible_pairs(M6, Gm, Gset):
    cp = cache_path(CACHE_INCOMPAT_FILE)
    if CACHE_PRECOMPUTE and os.path.exists(cp):
        log("Loading incompatibility from cache...")
        return load(cp)
    
    log(f"Computing incompatible pairs for {len(Gm)} flats...")
    t0 = time.time()
    
    pairs = [(i, j, int(Gm[i]) | int(Gm[j])) 
             for i in range(len(Gm)) for j in range(i+1, len(Gm))
             if not comparable(Gm[i], Gm[j])]
    
    log(f"  {len(pairs)} incomparable pairs to check")
    unions = list(set(p[2] for p in pairs))
    groundset = list(M6.groundset())
    closure_cache = {}
    
    for idx, mask in enumerate(unions):
        subset = [groundset[i] for i in mask_to_list(mask)]
        closure_cache[mask] = set_to_mask(M6.closure(subset))
        if (idx + 1) % 50000 == 0:
            log(f"    {idx+1}/{len(unions)} closures...")
    
    incompatible = {(i, j) for i, j, u in pairs if closure_cache[u] in Gset}
    log(f"  {len(incompatible)} incompatible pairs in {time.time()-t0:.1f}s")
    if CACHE_PRECOMPUTE:
        save(incompatible, cp)
    return incompatible

