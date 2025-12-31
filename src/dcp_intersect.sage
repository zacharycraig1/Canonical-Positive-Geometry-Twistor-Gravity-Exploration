# dcp_intersect.sage
# Assumes dcp_common.sage, dcp_search.sage, dcp_flats.sage are loaded

def _coeff_matrix_from_scan_B(B, scales):
    d = B.nrows()
    k = B.ncols()
    return matrix(QQ, d, k, lambda i,j: QQ(B[i,j]) * QQ(scales[i]))

def _combine_invariant_basis(Vinv_old, Cmat):
    d = len(Vinv_old)
    k = Cmat.ncols()
    zero = 0 * Vinv_old[0]
    Vinv_new = []
    for j in range(k):
        v = zero
        for i in range(d):
            cij = Cmat[i,j]
            if cij:
                v += cij * Vinv_old[i]
        Vinv_new.append(v)
    return Vinv_new

def _lr_to_json(lr):
    if lr is None: return None
    if isinstance(lr, (list, tuple)):
        if len(lr) == 0: return []
        first = lr[0]
        if isinstance(first, (list, tuple)) and len(first) == 2:
            out = []
            for a, b in lr: out.append([int(a), int(b)])
            return out
        out = []
        for x in lr:
            try: out.append(int(x))
            except Exception: out.append(x)
        return out
    try: return int(lr)
    except Exception: return lr

def intersection_across_boundaries(C, triples, triples_array, m, Gm, Gset, full_mask, incompatible, Vinv0, Vbasis, triple_sign, checkpoint_path=None):
    # Try to resume from checkpoint
    start_idx = 0
    if RESUME_FROM_CHECKPOINT and checkpoint_path and os.path.exists(checkpoint_path):
        try:
            log(f"Attempting to resume from checkpoint: {checkpoint_path}")
            checkpoint_data = load(checkpoint_path)
            start_idx = int(checkpoint_data.get('last_boundary_idx', 0))
            Vinv_cur = checkpoint_data.get('Vinv_cur')
            T = checkpoint_data.get('T')
            log_table = checkpoint_data.get('log_table', [])
            if Vinv_cur is not None and T is not None:
                Vinv_cur = list(Vinv_cur)
                log(f"  Resumed from boundary {start_idx}, current_dim={len(Vinv_cur)}")
            else:
                start_idx = 0
                Vinv_cur = list(Vinv0)
                d0 = len(Vinv0)
                T = identity_matrix(QQ, d0)
                log_table = []
        except Exception as e:
            log(f"  Checkpoint load failed ({e}), starting fresh")
            start_idx = 0
            Vinv_cur = list(Vinv0)
            d0 = len(Vinv0)
            T = identity_matrix(QQ, d0)
            log_table = []
    else:
        Vinv_cur = list(Vinv0)
        d0 = len(Vinv0)
        T = identity_matrix(QQ, d0)
        log_table = []

    for t, S in enumerate(INTERSECT_BOUNDARIES):
        if t < start_idx: continue
        d = len(Vinv_cur)
        log(f"\n{'='*70}\nINTERSECT {t+1}/{len(INTERSECT_BOUNDARIES)}  boundary S={S}  current_dim={d}\n{'='*70}")

        S_tuple = tuple(S)
        if len(S_tuple) == 2:
            S_ch = tuple(sorted(S_tuple))
            if S_ch not in C: continue
            Sidx = C.index(S_ch)
            Left = set(S_tuple)
            Right = set(range(1,7)) - Left
        elif len(S_tuple) == 3:
            S = canonical_three_subset(S_tuple, 6)
            S_ch = channel_of_three_subset(S, 6)
            Sidx = C.index(S_ch)
            Left = set(S)
            Right = set(range(1,7)) - Left
        else: continue

        Smask = PY1 << Sidx
        must_idx = [i for i, fm in enumerate(Gm) if int(fm) == Smask]
        Lonly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
        Ronly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)

        log(f"Boundary masks: Sidx={Sidx}, |Left|={len(Left)}, |Right|={len(Right)}")

        # Create a deterministic cache key for this boundary
        boundary_key = "_".join(map(str, sorted(list(S))))
        sols, best_size = run_search(Gm, Gset, full_mask, incompatible, must_idx, cache_key=boundary_key)
        if not sols:
            log("No charts found for this boundary. Skipping.")
            continue

        sol, sz = sols[0]
        res = scan_chart_exact_smallD(sol, Vinv_cur, triples, m, Lonly, Ronly, Smask, triple_sign, return_basis=True)
        status = res.get('status')
        null_dim = int(res.get('null_dim', 0))
        bad = res.get('bad', None)
        lr = res.get('lr', None)

        log(f"Result: status={status}  null_dim={null_dim}  bad={bad}  lr_total={lr}")

        log_table.append({
            'S': list(S), 'status': status, 'null_dim': null_dim,
            'bad': int(bad) if bad is not None else None,
            'lr': _lr_to_json(lr), 'chart_size': int(sz),
        })

        if status != 'HIT' or null_dim == 0:
            log("Intersection became empty on this boundary.")
            if INTERSECT_STOP_ON_EMPTY: break
            Vinv_cur = []
            T = matrix(QQ, T.nrows(), 0)
            break
        
        B = res.get('_B')
        scales = res.get('_scales')
        Cmat = _coeff_matrix_from_scan_B(B, scales)
        Vinv_cur = _combine_invariant_basis(Vinv_cur, Cmat)
        T = T * Cmat

        if checkpoint_path and (t + 1) % INTERSECT_CHECKPOINT_EVERY == 0:
            save({
                'last_boundary_idx': t + 1,
                'Vinv_cur': Vinv_cur,
                'T': T,
                'log_table': log_table,
            }, checkpoint_path)

        if len(Vinv_cur) == INTERSECT_TARGET_DIM:
            log(f"\n[INTERSECT] Reached target dimension {INTERSECT_TARGET_DIM}!")
            break

    alpha0 = None
    if len(Vinv_cur) > 0:
        alpha0 = T.column(0)
    
    return alpha0, log_table, T, Vinv_cur

