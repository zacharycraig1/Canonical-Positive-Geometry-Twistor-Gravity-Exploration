# dcp_common.sage
from sage.all import *
import numpy as np
import builtins as _builtins
import json
import gc
import multiprocessing as mp
import os, time, sys, shutil, itertools
import random
from collections import defaultdict
from itertools import combinations

# =============================================================================
# CONFIGURATION
# =============================================================================
DIAG = True

# Search settings
TOTAL_TRIALS = 50000
INTERSECTION_MODE = True
INTERSECT_BOUNDARIES = [
    (1,2,3), (1,2,4), (1,3,4), (2,3,4), (1,2,5), (1,4,5)
]
INTERSECT_BOUNDARY_MODE = "ALL_3x3"
INTERSECT_RANDOMIZE_ORDER = False
RUN_CANDIDATE_VALIDATION = True
VALIDATION_SHUFFLES = 1
VALIDATION_SEED_OFFSETS = [0]
VALIDATION_REQUIRE_DIM1 = True
SAVE_ARTIFACTS = True
ARTIFACT_DIR = "runs"
INTERSECT_MAX_CHARTS_PER_BOUNDARY = 1
INTERSECT_TARGET_DIM = 1
INTERSECT_STOP_ON_EMPTY = True
INTERSECT_CHECKPOINT_EVERY = 1
RESUME_FROM_CHECKPOINT = False

NUM_SEEDS = 1
BASE_SEED = 42
SEED_STRIDE = 1000
TRIALS_PER_SEED = TOTAL_TRIALS // NUM_SEEDS
MAX_SOLUTIONS = 100

CACHE_PRECOMPUTE = True
CACHE_CHARTS = True
CACHE_SCAN_PROGRESS = False

if CACHE_PRECOMPUTE:
    CACHE_DIR = os.path.join(os.getcwd(), 'dcp_cache')
else:
    CACHE_DIR = os.path.join(os.getcwd(), 'dcp_cache_' + time.strftime('%Y%m%d_%H%M%S'))

_TOTAL_CORES = os.cpu_count() or 2
SCAN_WORKERS = max(1, _TOTAL_CORES - 1)  # Leave one for OS

MODP_PREFILTER = True
MODP_SAMPLE_ROWS = 8000
MODP_PRIMES = (1000003, 1000033, 1000037, 1000039)

CACHE_FLATS_FILE = "flats.sobj"
CACHE_INCOMPAT_FILE = "incompat.sobj"
CACHE_SOLUTIONS_FILE = "solutions.sobj"
CACHE_SCAN_PROGRESS_FILE = "scan_progress.sobj"
CACHE_INVARIANTS_FILE = "invariants.sobj"

MAX_CHARTS_TO_SCAN = 10
CHECKPOINT_EVERY = 5
MP_CHUNKSIZE = 4
EXACT_SCAN = True
FORCE_FRESH = False
ONLY_BEST_SIZE = False
EARLY_ABORT_ON_CONTRADICTION = True
USE_PARI_KERNEL = True

INVARIANT_MODE = 'S3xS3'
MULTI_STRATEGY_SEARCH = True
STRATEGY_INVARIANT_MODES = ['S3xS3', 'S3xS3Z2', 'S6']
STRATEGY_STOP_ON_DIM1 = True

GLOBAL_VERIFY = True
GLOBAL_STOP_ON_FAIL = True
GLOBAL_MAX_BOUNDARIES = 20
WRITE_HIT_REPORT = True
HIT_REPORT_DIR = CACHE_DIR
PRINT_HIT_REPORT_TO_TERMINAL = True

PY0, PY1 = int(0), int(1)

# =============================================================================
# UTILITIES
# =============================================================================

def json_sanitize(x):
    try:
        from sage.all import Integer as SageInteger, Rational as SageRational
    except Exception:
        SageInteger = ()
        SageRational = ()
    if x is None or isinstance(x, (bool, int, float, str)):
        return x
    if SageInteger and isinstance(x, SageInteger):
        return int(x)
    if SageRational and isinstance(x, SageRational):
        return {'num': int(x.numerator()), 'den': int(x.denominator())}
    if isinstance(x, dict):
        return {str(k): json_sanitize(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [json_sanitize(v) for v in x]
    if isinstance(x, set):
        return [json_sanitize(v) for v in sorted(list(x))]
    try:
        return int(x)
    except Exception:
        return str(x)

def popcount(x): return int(x).bit_count()
def ts(): return time.strftime("%H:%M:%S")
def log(msg):
    if DIAG: print(f"[{ts()}] {msg}", flush=True)

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

def comparable(A, B):
    A, B = int(A), int(B)
    return ((A & ~B) == 0) or ((B & ~A) == 0)

def ensure_cache_dir():
    os.makedirs(CACHE_DIR, exist_ok=True)

def cache_path(name):
    ensure_cache_dir()
    return os.path.join(CACHE_DIR, name)

def clear_all_caches():
    if os.path.exists(CACHE_DIR):
        try:
            shutil.rmtree(CACHE_DIR)
            log(f"Cleared cache directory: {CACHE_DIR}")
        except Exception as e:
            log(f"Failed to clear cache: {e}")
    os.makedirs(CACHE_DIR, exist_ok=True)

def fast_right_kernel(M):
    if not USE_PARI_KERNEL:
        return M.right_kernel().basis_matrix()
    try:
        pari_result = M.__pari__().matker(1)
        if pari_result.ncols() == 0:
            return matrix(M.base_ring(), 0, M.ncols())
        kernel_basis = pari_result.mattranspose().sage()
        return matrix(M.base_ring(), kernel_basis)
    except Exception as e:
        log(f"    [PARI] fallback to standard kernel: {e}")
        return M.right_kernel().basis_matrix()

def fast_kernel_dim(M):
    if not USE_PARI_KERNEL:
        return M.right_kernel().dimension()
    try:
        return int(M.__pari__().matker(1).ncols())
    except:
        return M.right_kernel().dimension()

class SignTable:
    def __init__(self, max_size=600):
        log(f"Building sign lookup tables (max_size={max_size})...")
        t0 = time.time()
        self.max_size = max_size
        self.pair_sign = np.zeros((max_size, max_size), dtype=np.int8)
        for i in range(max_size):
            for j in range(i+1, max_size):
                self.pair_sign[i, j] = 1
                self.pair_sign[j, i] = -1
        log(f"  Sign tables built in {time.time()-t0:.2f}s")
    
    def get_pair_sign(self, i, j):
        if i < self.max_size and j < self.max_size:
            return int(self.pair_sign[i, j])
        if i < j: return 1
        if i > j: return -1
        return 0

SIGN_TABLE = None
def get_sign_table():
    global SIGN_TABLE
    if SIGN_TABLE is None:
        SIGN_TABLE = SignTable(600)
    return SIGN_TABLE

def triple_sign(u, v, w):
    """Parity sign of ordering of (u,v,w). Returns 0 if any repeats."""
    if u == v or u == w or v == w:
        return 0
    inv = 0
    if u > v: inv += 1
    if u > w: inv += 1
    if v > w: inv += 1
    return -1 if (inv & 1) else 1

def canonical_three_subset(S, n=6):
    S = set(S)
    comp = set(range(1, n+1)) - S
    return min(tuple(sorted(S)), tuple(sorted(comp)))

def channel_of_three_subset(S, n=6):
    if len(S) != 3:
        raise ValueError(f"channel_of_three_subset expects |S|=3")
    return canonical_three_subset(S, n)

def all_distinct_3x3_boundaries(n=6):
    out = []
    for S in itertools.combinations(range(1, n+1), 3):
        S = tuple(S)
        comp = tuple(sorted(set(range(1, n+1)) - set(S)))
        if S < comp:
            out.append(S)
    return out

def channels_C6():
    C = []
    for ij in Subsets(range(1,7), 2):
        C.append(tuple(sorted(ij)))
    for S in Subsets(range(1,7), 3):
        C.append(canonical_three_subset(S, 6))
    return sorted(set(C))

def ch_support(ch): return set(ch)

def build_M6_matroid():
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
    
    R = matrix(QQ, len(circuits4), Wdim, sparse=True)
    for r, cir in enumerate(circuits4):
        a, b, c, d = cir
        for sgn, tri in [(1,(b,c,d)), (-1,(a,c,d)), (1,(a,b,d)), (-1,(a,b,c))]:
            R[r, col_index[tuple(sorted(tri))]] += sgn
    
    pairs2 = [(i,j) for i in range(m) for j in range(i+1,m)]
    P = matrix(QQ, len(pairs2), Wdim, sparse=True)
    for r, (i,j) in enumerate(pairs2):
        for k in range(m):
            if k not in (i, j):
                sgn, tri = sign_sort((k, i, j))
                P[r, col_index[tri]] += sgn
    
    return triples, R.stack(P).right_kernel().basis_matrix()

def alpha_support(alpha):
    nz = []
    for i,a in enumerate(alpha):
        if a != 0:
            nz.append((i, a))
    return nz

def build_candidate_vec(alpha, Vinv):
    if alpha is None:
        return None
    Wdim = Vinv[0].nrows() if hasattr(Vinv[0],'nrows') else len(Vinv[0])
    v = vector(QQ, Wdim)
    for i,a in enumerate(alpha):
        if a:
            v += QQ(a) * vector(QQ, Vinv[i])
    return v

def _perm_ch(ch, p):
    if len(ch) == 2:
        a, b = p(ch[0]), p(ch[1])
        return (a, b) if a < b else (b, a)
    return canonical_three_subset({p(ch[0]), p(ch[1]), p(ch[2])}, 6)

def perm_for_three_subset(S):
    S = sorted(list(S))
    comp = sorted([x for x in range(1,7) if x not in S])
    imgs = S + comp
    return SymmetricGroup(6)(imgs)

def permute_flat_mask(Fmask, perm, C, chan_index):
    out = PY0
    m = int(Fmask)
    while m:
        lsb = m & -m
        idx = lsb.bit_length() - 1
        ch_new = _perm_ch(C[idx], perm)
        out |= (PY1 << chan_index[ch_new])
        m -= lsb
    return int(out)

def permute_chart(N, perm, C, chan_index):
    return [permute_flat_mask(F, perm, C, chan_index) for F in N]

def build_LR_masks(C, Left, Right):
    Lonly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Left), PY0)
    Ronly = _builtins.sum((PY1 << i for i, ch in enumerate(C) if ch_support(ch) <= Right), PY0)
    return int(Lonly), int(Ronly)

def solution_hash(sol):
    return tuple(sorted(int(x) for x in sol))
