
# src/dcp/proxy_pullback.sage
from sage.all import *

def scan_chart_proxy(N, Vinv, triples, m, Lonly, Ronly, Smask, triple_sign_func, return_details=False):
    """
    Proxy pullback scan (leading log approximation).
    Replicates scan_chart_exact2d but can return detailed residue info.
    """
    import time
    t0 = time.time()
    
    # 1. Setup
    N_sorted = sorted([int(x) for x in N], key=lambda x: (x.bit_count(), x))
    n_flats = len(N_sorted)
    N_index = {F: u for u, F in enumerate(N_sorted)}
    
    Smask_int = int(Smask)
    if Smask_int not in N_index:
        return {'status': 'skip', 'reason': 'no_S', 'time': time.time()-t0}
        
    uS = N_index[Smask_int]
    Sidx = Smask_int.bit_length() - 1
    
    # 2. Build phi_lists (which u-vars contribute to which channel)
    # phi_lists[ch] = list of u indices
    # Assumes u_F contributes to ch if ch is in F (and minimal?)
    # The proxy logic in dcp_search just checks "ch in mask_to_list(F)"
    # which implies u_F appears in ALL channels contained in F?
    # Wait, scan_chart_exact2d says:
    # for u, F in enumerate(N_sorted):
    #    for ch in mask_to_list(F): phi_lists[ch].append(u)
    # mask_to_list(F) returns list of elements (channels) in F.
    # So yes, u_F appears in every channel variable that is in F.
    
    phi_lists = [[] for _ in range(m)]
    
    def mask_to_list(mask):
        out, mm = [], int(mask)
        while mm:
            lsb = mm & -mm
            out.append(lsb.bit_length() - 1)
            mm -= lsb
        return out
        
    for u, F in enumerate(N_sorted):
        for ch in mask_to_list(F):
            phi_lists[ch].append(u)
            
    # 3. L/R Classification of u variables
    cls = [0] * n_flats
    Lonly_int = int(Lonly)
    Ronly_int = int(Ronly)
    
    for u, F in enumerate(N_sorted):
        if F != 0:
            if (F & ~Lonly_int) == 0:
                cls[u] = 1 # Left
            elif (F & ~Ronly_int) == 0:
                cls[u] = 2 # Right
                
    # 4. Process Invariants
    # Vinv is list of vectors (dictionaries)
    d = len(Vinv)
    vinv_data = []
    scales = []
    
    for vec in Vinv:
        # If vec is a vector object, get dict, else assume dict
        items = vec.dict().items() if hasattr(vec, 'dict') else vec.items()
        kept_raw = []
        denoms = []
        for k, v in items:
            k = int(k)
            # Filter terms that don't involve S (uS)
            # The residue calculation requires at least one channel in the triple to produce uS pole?
            # scan_chart_exact2d logic:
            ci, cj, ck = triples[k]
            if (ci != Sidx) and (cj != Sidx) and (ck != Sidx):
                continue
                
            vq = QQ(v)
            kept_raw.append((k, vq))
            denoms.append(int(vq.denominator()))
            
        if not denoms:
            scale = 1
        else:
            scale = 1
            for dd in denoms: scale = lcm(scale, dd)
            
        scales.append(QQ(scale))
        kept = []
        for k, vq in kept_raw:
            vi = vq * scale
            if vi.denominator() == 1:
                kept.append((k, int(vi)))
            else:
                kept.append((k, QQ(vi)))
        vinv_data.append(kept)
        
    # 5. Compute Residues
    # We want to find linear combination alpha such that Residue is "good" (L-only * R-only)
    # Good terms: (L, R) or (R, L)
    # Bad terms: (L, L), (R, R), (M, X), (X, M), (X, X) where X is mixed/neither
    # Wait, scan_chart_exact2d counts "lr_count" for (1,2) or (2,1).
    # Everything else is "bad" and put into bad_coeffs.
    
    from collections import defaultdict
    bad_coeffs = defaultdict(lambda: [0]*d)
    lr_count = [0]*d
    
    # We process each basis vector separately
    for v_idx in range(d):
        for col_idx, coeff in vinv_data[v_idx]:
            ci, cj, ck = triples[int(col_idx)]
            pi = phi_lists[ci]
            pj = phi_lists[cj]
            pk = phi_lists[ck]
            
            if not pi or not pj or not pk: continue
            
            # Logic from scan_chart_exact2d
            # We are looking for poles at uS=0.
            # dlog(l_ci) ^ dlog(l_cj) ^ dlog(l_ck)
            # Each l_x is product of u_y.
            # Res_{uS=0} is sum over terms where exactly one of l_ci, l_cj, l_ck contains uS,
            # and we drop that dlog uS and keep the others.
            
            # Cases:
            # 1. ci has uS, others don't?
            # Actually, scan_chart_exact2d iterates over all u,v,w in pi,pj,pk
            # and checks if one of them is uS.
            
            # Optimization: if ci==Sidx, then l_ci MUST contain uS (since S is in Chain(S))
            # But generally, uS is in l_x if S <= flat(x) ?? No.
            # In proxy: l_e = product u_F where e in F.
            # So uS is in l_e if e in S.
            # S is a flat. mask_to_list(S) gives channels in S.
            # So if e in S, uS is in l_e.
            
            # scan_chart_exact2d logic:
            
            if ci == Sidx:
                # uS is in l_ci.
                # residue is dlog(l_cj)^dlog(l_ck) * sgn
                u = uS
                for v in pj:
                    for w in pk:
                        # Term dlog u ^ dlog v ^ dlog w.
                        # Res_u=0 -> dlog v ^ dlog w.
                        a, b = v, w
                        ca, cb = cls[a], cls[b]
                        if (ca == 1 and cb == 1) or (ca == 2 and cb == 2): continue
                        
                        sgn = triple_sign_func(u, v, w)
                        if sgn:
                            # normalize pair (a,b)
                            if a > b: pair = (b, a); sgn = -sgn
                            else: pair = (a, b)
                            
                            is_lr = (cls[pair[0]] == 1 and cls[pair[1]] == 2) or \
                                    (cls[pair[0]] == 2 and cls[pair[1]] == 1)
                                    
                            if is_lr:
                                lr_count[v_idx] += 1
                            else:
                                bad_coeffs[pair][v_idx] += sgn * coeff
                                
            elif cj == Sidx:
                v = uS
                for u in pi:
                    for w in pk:
                        a, b = u, w
                        ca, cb = cls[a], cls[b]
                        if (ca == 1 and cb == 1) or (ca == 2 and cb == 2): continue
                        sgn = triple_sign_func(u, v, w)
                        if sgn:
                            if a > b: pair = (b, a); sgn = -sgn
                            else: pair = (a, b)
                            is_lr = (cls[pair[0]] == 1 and cls[pair[1]] == 2) or \
                                    (cls[pair[0]] == 2 and cls[pair[1]] == 1)
                            if is_lr: lr_count[v_idx] += 1
                            else: bad_coeffs[pair][v_idx] += sgn * coeff

            elif ck == Sidx:
                w = uS
                for u in pi:
                    for v in pj:
                        a, b = u, v
                        ca, cb = cls[a], cls[b]
                        if (ca == 1 and cb == 1) or (ca == 2 and cb == 2): continue
                        sgn = triple_sign_func(u, v, w)
                        if sgn:
                            if a > b: pair = (b, a); sgn = -sgn
                            else: pair = (a, b)
                            is_lr = (cls[pair[0]] == 1 and cls[pair[1]] == 2) or \
                                    (cls[pair[0]] == 2 and cls[pair[1]] == 1)
                            if is_lr: lr_count[v_idx] += 1
                            else: bad_coeffs[pair][v_idx] += sgn * coeff
            
            else:
                # General case: iterate all, check if uS is present
                for u in pi:
                    for v in pj:
                        for w in pk:
                            if u == uS:
                                a, b = v, w
                            elif v == uS:
                                a, b = u, w
                            elif w == uS:
                                a, b = u, v
                            else:
                                continue
                                
                            ca, cb = cls[a], cls[b]
                            if (ca == 1 and cb == 1) or (ca == 2 and cb == 2): continue
                            sgn = triple_sign_func(u, v, w)
                            if sgn:
                                if a > b: pair = (b, a); sgn = -sgn
                                else: pair = (a, b)
                                is_lr = (cls[pair[0]] == 1 and cls[pair[1]] == 2) or \
                                        (cls[pair[0]] == 2 and cls[pair[1]] == 1)
                                if is_lr: lr_count[v_idx] += 1
                                else: bad_coeffs[pair][v_idx] += sgn * coeff

    # 6. Analyze Results
    bad_list = [(pair, val) for pair, val in bad_coeffs.items() if any(v != 0 for v in val)]
    n_bad = len(bad_list)
    
    result = {
        'n_flats': n_flats,
        'bad_count': n_bad,
        'lr_count': lr_count,
        'scales': scales
    }
    
    if return_details:
        result['bad_coeffs'] = bad_list
        result['N'] = N
        result['u_map'] = {u: F for u, F in enumerate(N_sorted)}
        result['cls'] = cls
    
    # Solve for alpha
    if d == 2:
        # Special case for dim 2
        # We want alpha[0]*bad[pair][0] + alpha[1]*bad[pair][1] = 0 for all pairs
        
        # Simple solver for 2 vars
        sol = None
        if n_bad == 0:
            sol = (1, 0) # Arbitrary
        else:
            # Pick first non-zero constraint
            for pair, coeffs in bad_list:
                c0, c1 = coeffs
                if c0 == 0 and c1 == 0: continue
                if c0 == 0:
                    # c1 * a1 = 0 => a1=0. a0=1
                    sol = (1, 0)
                elif c1 == 0:
                    # c0 * a0 = 0 => a0=0. a1=1
                    sol = (0, 1)
                else:
                    # c0*a0 + c1*a1 = 0 => a0/a1 = -c1/c0
                    # a0 = c1, a1 = -c0
                    sol = (c1, -c0)
                break
                
        # Verify solution against all other constraints
        valid = True
        contradictions = 0
        if sol:
             a0, a1 = sol
             for pair, coeffs in bad_list:
                 if coeffs[0]*a0 + coeffs[1]*a1 != 0:
                     valid = False
                     contradictions += 1
                     
        result['alpha_sol'] = sol
        result['valid'] = valid
        result['contradictions'] = contradictions
        
        if valid:
            result['status'] = 'HIT'
        else:
            result['status'] = 'NO_HIT'
            
    else:
        # General solver (not implemented here for >2, using placeholder)
        result['status'] = 'UNKNOWN_DIM'
        
    return result







