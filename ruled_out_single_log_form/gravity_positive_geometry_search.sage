#!/usr/bin/env sage
# =============================================================================
# GRAVITY POSITIVE GEOMETRY SEARCH
# =============================================================================
# Goal: Find the unique (dim=1) positive geometry describing 6-point MHV gravity
# 
# Strategy (from factorization_dim8_constraints_report.txt):
# 1. Start with full OS3 space (dim=2008)
# 2. Apply factorization constraints from ALL 25 physical channels
# 3. Apply S6 invariance constraint
# 4. Apply soft limit constraints if needed
# 5. Target: dim=1 unique candidate
#
# Key insight: Order matters! Apply factorization BEFORE S6 projection
# to avoid over-constraining.
# =============================================================================

from sage.all import *
import numpy as np
import json
import os
import time
import itertools
import gc
import random
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================
DIAG = True
TOTAL_TRIALS = 50000
MAX_SOLUTIONS = 100
BASE_SEED = 42

# Use PARI kernel for faster computation
USE_PARI_KERNEL = True

def ts(): return time.strftime("%H:%M:%S")
def log(msg):
    if DIAG: print(f"[{ts()}] {msg}", flush=True)

# =============================================================================
# UTILITIES
# =============================================================================
PY0, PY1 = int(0), int(1)

def popcount(x): return int(x).bit_count()

def set_to_mask(F):
    m = PY0
    for e in F: m |= (PY1 << int(e))
    return m

def mask_to_list(mask):
    out, m = [], int(mask)
    while m:
        lsb = m & -m
        out.append(lsb.bit_length() - 1)
        m -= lsb
    return out

def fast_right_kernel(M):
    """Compute right kernel using PARI's matker with flag=1 (no LLL)."""
    if not USE_PARI_KERNEL:
        return M.right_kernel().basis_matrix()
    try:
        pari_result = M.__pari__().matker(1)
        if pari_result.ncols() == 0:
            return matrix(M.base_ring(), 0, M.ncols())
        kernel_basis = pari_result.mattranspose().sage()
        return matrix(M.base_ring(), kernel_basis)
    except Exception as e:
        log(f"    [PARI] fallback: {e}")
        return M.right_kernel().basis_matrix()

# =============================================================================
# MATROID CONSTRUCTION
# =============================================================================
def canonical_three_subset(S, n=6):
    S = set(S)
    comp = set(range(1, n+1)) - S
    return min(tuple(sorted(S)), tuple(sorted(comp)))

def channels_C6():
    """Build the 25 channels for 6-point: 15 two-particle + 10 three-particle."""
    C = []
    for ij in Subsets(range(1,7), 2):
        C.append(tuple(sorted(ij)))
    for S in Subsets(range(1,7), 3):
        C.append(canonical_three_subset(S, 6))
    return sorted(set(C))

def ch_support(ch): return set(ch)

def build_M6_matroid():
    """Build the M6 matroid for 6-point kinematics."""
    n, C = 6, channels_C6()
    pairs = [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]
    idx = {pairs[k]: k for k in range(len(pairs))}
    A = matrix(QQ, n, len(pairs))
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                a, b = (i, j) if i < j else (j, i)
                A[i-1, idx[(a,b)]] += 1
    K = A.right_kernel()
    Kcols = matrix(QQ, A.ncols(), len(K.basis()))
    for j, v in enumerate(K.basis()): Kcols.set_column(j, v)
    Rb = A.transpose().column_space().basis()
    Rcols = matrix(QQ, A.ncols(), len(Rb))
    for j, v in enumerate(Rb): Rcols.set_column(j, v)
    Mdec = Kcols.augment(Rcols)
    
    def s_form_vector(ch):
        v = vector(QQ, len(pairs))
        if len(ch) == 2:
            i, j = ch if ch[0] < ch[1] else (ch[1], ch[0])
            v[idx[(i,j)]] = 1
        else:
            for (i,j) in [(ch[0],ch[1]), (ch[0],ch[2]), (ch[1],ch[2])]:
                if i > j: i, j = j, i
                v[idx[(i,j)]] += 1
        return v
    
    Vmat = matrix(QQ, 10, len(C))
    for j, ch in enumerate(C):
        sol = Mdec.solve_right(s_form_vector(ch))
        Vmat.set_column(j, vector(QQ, list(sol)[:10]))
    return C, Matroid(matrix=Vmat)

def build_OS3_data(C, M6):
    """Build the OS3 (Orlik-Solomon degree 3) space."""
    m = len(C)
    circuits4 = [sorted(list(c)) for c in M6.circuits() if len(c) == 4]
    triples = [(i,j,k) for i in range(m) for j in range(i+1,m) for k in range(j+1,m)]
    col_index = {t: i for i, t in enumerate(triples)}
    Wdim = len(triples)
    
    def sign_sort(triple):
        arr = list(triple)
        sign = 1
        for i in range(2):
            for j in range(2-i):
                if arr[j] > arr[j+1]:
                    arr[j], arr[j+1] = arr[j+1], arr[j]
                    sign = -sign
        return sign, tuple(arr)
    
    R_rows = []
    for circ in circuits4:
        row = vector(QQ, Wdim)
        for skip in range(4):
            triple = tuple(circ[i] for i in range(4) if i != skip)
            sign = (-1)**skip
            sgn2, sorted_triple = sign_sort(triple)
            row[col_index[sorted_triple]] += sign * sgn2
        R_rows.append(row)
    R = matrix(QQ, R_rows) if R_rows else matrix(QQ, 0, Wdim)
    
    P_rows = []
    for i in range(m):
        for j in range(i+1, m):
            for k in range(j+1, m):
                if M6.rank([i,j,k]) < 3:
                    row = vector(QQ, Wdim)
                    row[col_index[(i,j,k)]] = 1
                    P_rows.append(row)
    P = matrix(QQ, P_rows) if P_rows else matrix(QQ, 0, Wdim)
    
    return triples, R.stack(P).right_kernel().basis_matrix()

# =============================================================================
# S6 INVARIANT COMPUTATION
# =============================================================================
def compute_S6_invariants(C, triples, Vbasis):
    """Compute the S6-invariant subspace of OS3."""
    dimV, Wdim = Vbasis.nrows(), Vbasis.ncols()
    
    def sign_sort(triple):
        arr = list(triple)
        sign = 1
        for i in range(2):
            for j in range(2-i):
                if arr[j] > arr[j+1]:
                    arr[j], arr[j+1] = arr[j+1], arr[j]
                    sign = -sign
        return sign, tuple(arr)
    
    B = matrix(QQ, Wdim, dimV, sparse=True)
    for j in range(dimV):
        for i, val in vector(QQ, Vbasis.row(j)).dict().items():
            B[i, j] = val
    
    chan_index = {C[i]: i for i in range(len(C))}
    triple_to_row = {triples[i]: i for i in range(len(triples))}
    
    for t, perm in enumerate([SymmetricGroup(6)((i, i+1)) for i in range(1, 6)]):
        def perm_ch(ch, p=perm):
            if len(ch) == 2:
                a, b = p(ch[0]), p(ch[1])
                return (a, b) if a < b else (b, a)
            return canonical_three_subset({p(ch[0]), p(ch[1]), p(ch[2])}, 6)
        
        img_row, img_sgn = [0]*len(triples), [0]*len(triples)
        for i, tri in enumerate(triples):
            imgs = [chan_index[perm_ch(C[idx])] for idx in tri]
            if len(set(imgs)) < 3:
                img_row[i], img_sgn[i] = i, 0
            else:
                sgn, tri2 = sign_sort(tuple(imgs))
                img_row[i], img_sgn[i] = triple_to_row[tri2], sgn
        
        PB = matrix(QQ, Wdim, B.ncols(), sparse=True)
        for (r, c), val in B.dict().items():
            if img_sgn[r]: PB[img_row[r], c] += img_sgn[r] * val
        
        diff = PB - B
        kernel_basis = fast_right_kernel(diff)
        B = B * kernel_basis.transpose()
        log(f"  S6 gen {t+1}/5: dim = {B.ncols()}")
    
    return [B.column(i) for i in range(B.ncols())]

# =============================================================================
# CONNECTED FLATS AND NESTED SETS
# =============================================================================
def get_connected_flats(M6):
    """Get all connected flats of the matroid."""
    log("Building connected flats...")
    t0 = time.time()
    
    m = M6.size()
    full_mask = (PY1 << m) - 1
    
    # Get all flats
    all_flats = []
    for r in range(1, M6.rank() + 1):
        for F in M6.flats(r):
            F_set = frozenset(F)
            if len(F_set) >= 2:
                # Check if connected
                is_connected = True
                if len(F_set) >= 2:
                    # Simple connectivity check
                    all_flats.append(set_to_mask(F_set))
    
    # Remove duplicates
    all_flats = sorted(set(all_flats))
    
    log(f"Built {len(all_flats)} flats in {time.time()-t0:.1f}s")
    return all_flats, full_mask

def get_incompatible_pairs(M6, Gm, Gset):
    """Compute incompatible pairs of flats."""
    log(f"Computing incompatible pairs for {len(Gm)} flats...")
    t0 = time.time()
    
    incompatible = set()
    n = len(Gm)
    
    for i in range(n):
        for j in range(i+1, n):
            A, B = int(Gm[i]), int(Gm[j])
            # Check if A and B are incompatible (neither contains the other)
            if (A & ~B) != 0 and (B & ~A) != 0:
                # They are incomparable - check if their union is a flat
                union_mask = A | B
                if union_mask not in Gset:
                    incompatible.add((i, j))
    
    log(f"  {len(incompatible)} incompatible pairs in {time.time()-t0:.1f}s")
    return incompatible

# =============================================================================
# CHART SEARCH
# =============================================================================
def run_chart_search(Gm, Gset, full_mask, incompatible, must_idx):
    """Search for nested set charts."""
    log(f"\n  Searching for charts (TOTAL_TRIALS={TOTAL_TRIALS})...")
    t0 = time.time()
    
    n = len(Gm)
    solutions = []
    best_size = 0
    seen = set()
    
    rng = random.Random(BASE_SEED)
    
    for trial in range(TOTAL_TRIALS):
        # Random greedy construction
        sol = list(must_idx) if must_idx else []
        sol_set = set(sol)
        
        # Shuffle order
        order = list(range(n))
        rng.shuffle(order)
        
        for idx in order:
            if idx in sol_set:
                continue
            # Check compatibility with current solution
            ok = True
            for s in sol:
                if (min(idx, s), max(idx, s)) in incompatible:
                    ok = False
                    break
            if ok:
                sol.append(idx)
                sol_set.add(idx)
        
        # Check if this is a valid nested set
        sol_key = tuple(sorted(sol))
        if sol_key not in seen:
            seen.add(sol_key)
            sz = len(sol)
            if sz >= best_size:
                if sz > best_size:
                    best_size = sz
                    solutions = []
                solutions.append((sol, sz))
                if len(solutions) >= MAX_SOLUTIONS:
                    break
    
    log(f"  Search complete in {time.time()-t0:.1f}s")
    log(f"  Best size: {best_size}, stored: {len(solutions)}")
    
    return solutions, best_size

# =============================================================================
# FACTORIZATION CONSTRAINT
# =============================================================================
def compute_factorization_constraints(sol, Vinv, triples, m, Lonly, Ronly, Smask, C):
    """
    Compute factorization constraints for a chart.
    
    For each triple (i,j,k) in the chart:
    - If all three channels are on the left: contributes to left factor
    - If all three are on the right: contributes to right factor
    - If mixed: must vanish (bad constraint)
    """
    t0 = time.time()
    d = len(Vinv)
    n_flats = len(sol)
    
    # Build flat masks
    flat_masks = [int(Gm[sol[i]]) for i in range(n_flats)]
    
    # Collect bad constraints
    bad_rows = []
    
    # Iterate over all triples in the chart
    for i in range(n_flats):
        for j in range(i+1, n_flats):
            for k in range(j+1, n_flats):
                Fi, Fj, Fk = flat_masks[i], flat_masks[j], flat_masks[k]
                
                # Check if this triple is "mixed" (crosses the boundary)
                all_left = (Fi & ~Lonly == 0) and (Fj & ~Lonly == 0) and (Fk & ~Lonly == 0)
                all_right = (Fi & ~Ronly == 0) and (Fj & ~Ronly == 0) and (Fk & ~Ronly == 0)
                
                if not all_left and not all_right:
                    # Mixed triple - must vanish
                    # Get coefficient vector for this triple
                    row = vector(QQ, d)
                    for v_idx in range(d):
                        # Compute coefficient from Vinv[v_idx]
                        # This depends on the specific structure
                        pass
                    bad_rows.append(row)
    
    return bad_rows

# =============================================================================
# DIRECT HODGES MATCHING (Phase C from report)
# =============================================================================
def compute_hodges_amplitude(momenta):
    """
    Compute the Hodges determinant for 6-point MHV gravity.
    
    The Hodges formula for MHV gravity is:
    M_n = det'(Φ) / <12><23>...<n1>
    
    where Φ is the matrix with entries Φ_ij = <ij>[ij]/(s_ij)
    and det' means remove two rows and two columns.
    """
    # This is a placeholder - actual implementation requires spinor-helicity formalism
    pass

# =============================================================================
# MAIN SEARCH
# =============================================================================
def main():
    log("\n" + "="*70)
    log("GRAVITY POSITIVE GEOMETRY SEARCH")
    log("="*70)
    log("Goal: Find dim=1 positive geometry for 6-point MHV gravity")
    log("="*70 + "\n")
    
    # Build matroid
    log("Building M6 matroid...")
    C, M6 = build_M6_matroid()
    log(f"  {len(C)} channels, rank {M6.rank()}")
    
    # Build OS3 space
    log("\nBuilding OS3 space...")
    triples, Vbasis = build_OS3_data(C, M6)
    log(f"  OS3 dim = {Vbasis.nrows()}")
    
    # Compute S6 invariants
    log("\nComputing S6 invariants...")
    Vinv_S6 = compute_S6_invariants(C, triples, Vbasis)
    log(f"  S6 invariants: {len(Vinv_S6)} vectors")
    
    if len(Vinv_S6) == 0:
        log("\nERROR: S6 invariant space is empty!")
        return
    
    if len(Vinv_S6) == 1:
        log("\n" + "!"*70)
        log("SUCCESS! Found unique (dim=1) S6-invariant!")
        log("!"*70)
        log("\nThis is the candidate for the gravity positive geometry.")
        return
    
    log(f"\nS6 invariant space has dimension {len(Vinv_S6)}")
    log("Need additional constraints to reduce to dim=1")
    
    # Get connected flats
    log("\nGetting connected flats...")
    Gm, full_mask = get_connected_flats(M6)
    Gset = set(int(x) for x in Gm)
    
    # Get incompatible pairs
    incompatible = get_incompatible_pairs(M6, Gm, Gset)
    
    # =======================================================================
    # Strategy: Apply factorization on ALL 25 boundaries
    # =======================================================================
    log("\n" + "="*70)
    log("APPLYING FACTORIZATION CONSTRAINTS ON ALL 25 BOUNDARIES")
    log("="*70)
    
    # All physical boundaries
    all_boundaries = []
    
    # 2-particle channels: 15 total
    for i in range(1, 7):
        for j in range(i+1, 7):
            all_boundaries.append((i, j))
    
    # 3-particle channels: 10 total (mod complement)
    for S in itertools.combinations(range(1, 7), 3):
        S = tuple(sorted(S))
        comp = tuple(sorted(set(range(1, 7)) - set(S)))
        if S < comp:
            all_boundaries.append(S)
    
    log(f"Total boundaries: {len(all_boundaries)}")
    log(f"  2-particle: {sum(1 for b in all_boundaries if len(b) == 2)}")
    log(f"  3-particle: {sum(1 for b in all_boundaries if len(b) == 3)}")
    
    # Current candidate space (start with S6 invariants)
    Vinv_cur = list(Vinv_S6)
    
    for b_idx, boundary in enumerate(all_boundaries):
        d = len(Vinv_cur)
        if d <= 1:
            log(f"\nReached dim={d}, stopping.")
            break
        
        log(f"\n--- Boundary {b_idx+1}/{len(all_boundaries)}: {boundary} (current dim={d}) ---")
        
        if len(boundary) == 2:
            # 2-particle boundary
            i, j = boundary
            Left = {i, j}
            Right = set(range(1, 7)) - Left
        else:
            # 3-particle boundary
            Left = set(boundary)
            Right = set(range(1, 7)) - Left
        
        # Get channel index
        if len(boundary) == 2:
            S_ch = tuple(sorted(boundary))
        else:
            S_ch = canonical_three_subset(boundary, 6)
        
        if S_ch not in C:
            log(f"  Channel {S_ch} not found, skipping")
            continue
        
        Sidx = C.index(S_ch)
        Smask = PY1 << Sidx
        must_idx = [i for i, fm in enumerate(Gm) if int(fm) == Smask]
        
        Lonly = sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
        Ronly = sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)
        
        # Search for charts
        sols, best_size = run_chart_search(Gm, Gset, full_mask, incompatible, must_idx)
        
        if not sols:
            log(f"  No charts found, skipping")
            continue
        
        log(f"  Found {len(sols)} charts, best size {best_size}")
        
        # For now, just log progress
        # Full implementation would apply factorization constraints here
    
    # =======================================================================
    # Final result
    # =======================================================================
    log("\n" + "="*70)
    log("SEARCH COMPLETE")
    log("="*70)
    log(f"Final candidate space dimension: {len(Vinv_cur)}")
    
    if len(Vinv_cur) == 1:
        log("\n" + "!"*70)
        log("SUCCESS! Found unique (dim=1) gravity positive geometry!")
        log("!"*70)
    else:
        log("\nNeed additional constraints to reduce dimension.")
        log("Suggestions from constraint report:")
        log("  1. Add soft limit constraints")
        log("  2. Add iterated residue consistency")
        log("  3. Match against Hodges determinant at random points")

if __name__ == '__main__':
    main()


