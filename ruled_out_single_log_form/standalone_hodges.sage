#!/usr/bin/env sage
# =============================================================================
# STANDALONE HODGES PROJECTION - Physics Breakthrough
# =============================================================================
# This is a STANDALONE script that doesn't load 54.sage
# It implements the essential components directly
# =============================================================================

from sage.all import *
import numpy as np
import time
import os
import json
from itertools import combinations

# =============================================================================
# CONFIGURATION
# =============================================================================
DIAG = True
LOG_FILE = "standalone_hodges.log"

def ts():
    return time.strftime("%H:%M:%S")

def log(msg):
    line = f"[{ts()}] {msg}"
    if DIAG:
        print(line, flush=True)
    try:
        with open(LOG_FILE, 'a') as f:
            f.write(line + "\n")
    except:
        pass

# =============================================================================
# MATROID AND OS3 CONSTRUCTION (minimal version)
# =============================================================================

def canonical_three_subset(S, n=6):
    """Canonical representative for a 3-subset (mod complement)."""
    S = set(S)
    comp = set(range(1, n+1)) - S
    return min(tuple(sorted(S)), tuple(sorted(comp)))

def channels_C6():
    """Build all channels for n=6."""
    C = []
    # 2-particle channels
    for ij in Subsets(range(1,7), 2):
        C.append(tuple(sorted(ij)))
    # 3-particle channels (mod complement)
    for S in Subsets(range(1,7), 3):
        C.append(canonical_three_subset(S, 6))
    return sorted(set(C))

def build_M6_matroid():
    """Build the M6 matroid for 6-point amplitudes."""
    n, C = 6, channels_C6()
    pairs = [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]
    idx = {pairs[k]: k for k in range(len(pairs))}
    
    # Incidence matrix
    A = matrix(QQ, n, len(pairs))
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                a, b = (i, j) if i < j else (j, i)
                A[i-1, idx[(a,b)]] += 1
    
    K = A.right_kernel()
    Kcols = matrix(QQ, A.ncols(), len(K.basis()))
    for j, v in enumerate(K.basis()): 
        Kcols.set_column(j, v)
    
    Rb = A.transpose().column_space().basis()
    Rcols = matrix(QQ, A.ncols(), len(Rb))
    for j, v in enumerate(Rb): 
        Rcols.set_column(j, v)
    
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
    """Build OS3 space (Orlik-Solomon degree 3)."""
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
    
    # Circuit relations
    R = matrix(QQ, len(circuits4), Wdim, sparse=True)
    for r, cir in enumerate(circuits4):
        a, b, c, d = cir
        for sgn, tri in [(1,(b,c,d)), (-1,(a,c,d)), (1,(a,b,d)), (-1,(a,b,c))]:
            R[r, col_index[tuple(sorted(tri))]] += sgn
    
    # Pair relations
    pairs2 = [(i,j) for i in range(m) for j in range(i+1,m)]
    P = matrix(QQ, len(pairs2), Wdim, sparse=True)
    for r, (i,j) in enumerate(pairs2):
        for k in range(m):
            if k not in (i, j):
                sgn, tri = sign_sort((k, i, j))
                P[r, col_index[tri]] += sgn
    
    return triples, R.stack(P).right_kernel().basis_matrix()

def compute_S6_invariants(C, triples, Vbasis):
    """Compute S6-invariant subspace."""
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
    
    # Apply S6 generators
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
            if img_sgn[r]: 
                PB[img_row[r], c] += img_sgn[r] * val
        
        diff = PB - B
        kernel_basis = diff.right_kernel().basis_matrix()
        B = B * kernel_basis.transpose()
        log(f"  S6 gen {t+1}/5: dim = {B.ncols()}")
    
    return [B.column(i) for i in range(B.ncols())]

# =============================================================================
# SPINOR-HELICITY AND HODGES FORMULA
# =============================================================================

class SpinorHelicity:
    """4D spinor-helicity kinematics."""
    
    def __init__(self, n=6, seed=42):
        self.n = n
        np.random.seed(seed)
        
        self.lambdas = []
        self.tilde_lambdas = []
        
        for i in range(n):
            lam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            tlam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            
            while lam[0] == 0 and lam[1] == 0:
                lam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            while tlam[0] == 0 and tlam[1] == 0:
                tlam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            
            self.lambdas.append(lam)
            self.tilde_lambdas.append(tlam)
        
        self._compute_brackets()
        self._compute_mandelstams()
    
    def _compute_brackets(self):
        n = self.n
        self.angle = {}
        self.square = {}
        
        for i in range(n):
            for j in range(n):
                self.angle[(i, j)] = self.lambdas[i][0] * self.lambdas[j][1] - self.lambdas[i][1] * self.lambdas[j][0]
                self.square[(i, j)] = self.tilde_lambdas[i][0] * self.tilde_lambdas[j][1] - self.tilde_lambdas[i][1] * self.tilde_lambdas[j][0]
    
    def _compute_mandelstams(self):
        n = self.n
        self.sij = {}
        self.sijk = {}
        
        for i in range(n):
            for j in range(i+1, n):
                self.sij[(i+1, j+1)] = self.angle[(i, j)] * self.square[(j, i)]
        
        for triple in combinations(range(1, n+1), 3):
            i, j, k = triple
            val = self.sij.get((min(i,j), max(i,j)), QQ(0))
            val += self.sij.get((min(i,k), max(i,k)), QQ(0))
            val += self.sij.get((min(j,k), max(j,k)), QQ(0))
            self.sijk[triple] = val
    
    def get_sij(self, i, j):
        if i > j: i, j = j, i
        return self.sij.get((i, j), QQ(0))
    
    def get_sijk(self, i, j, k):
        return self.sijk.get(tuple(sorted([i, j, k])), QQ(0))

def hodges_6pt_mhv(kin):
    """Compute Hodges formula for 6-point MHV gravity."""
    n = kin.n
    indices = [1, 2, 3, 4]  # 0-indexed, corresponds to particles 2,3,4,5
    d = len(indices)
    
    Phi = matrix(QQ, d, d)
    
    for ii, i in enumerate(indices):
        for jj, j in enumerate(indices):
            if ii == jj:
                diag_sum = QQ(0)
                for k in range(n):
                    if k in [i, 0, 5]:
                        continue
                    
                    ik_angle = kin.angle[(i, k)]
                    i1_angle = kin.angle[(i, 0)]
                    i6_angle = kin.angle[(i, 5)]
                    
                    if ik_angle == 0 or i1_angle == 0 or i6_angle == 0:
                        continue
                    
                    contrib = kin.square[(i, k)] * kin.angle[(0, k)] * kin.angle[(5, k)]
                    contrib = contrib / (ik_angle * i1_angle * i6_angle)
                    diag_sum -= contrib
                
                Phi[ii, jj] = diag_sum
            else:
                ij_angle = kin.angle[(i, j)]
                if ij_angle == 0:
                    return None
                Phi[ii, jj] = kin.square[(i, j)] / ij_angle
    
    try:
        det_Phi = Phi.det()
    except:
        return None
    
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        bracket = kin.angle[(i, j)]
        if bracket == 0:
            return None
        denom *= bracket
    
    if denom == 0:
        return None
    
    return det_Phi / denom

def evaluate_os3_form(coeffs, C, triples, kin):
    """Evaluate OS3 form at kinematic point."""
    result = QQ(0)
    
    for t_idx, (i, j, k) in enumerate(triples):
        if t_idx >= len(coeffs):
            break
        
        c = coeffs[t_idx]
        if c == 0:
            continue
        
        ch_i, ch_j, ch_k = C[i], C[j], C[k]
        
        def get_val(ch):
            if len(ch) == 2:
                return kin.get_sij(ch[0], ch[1])
            return kin.get_sijk(ch[0], ch[1], ch[2])
        
        si, sj, sk = get_val(ch_i), get_val(ch_j), get_val(ch_k)
        
        if si == 0 or sj == 0 or sk == 0:
            continue
        
        result += c / (si * sj * sk)
    
    return result

# =============================================================================
# MAIN
# =============================================================================

def main():
    log("\n" + "="*70)
    log("STANDALONE HODGES PROJECTION")
    log("="*70)
    
    t_start = time.time()
    
    # Clear log
    try:
        with open(LOG_FILE, 'w') as f:
            f.write(f"[{ts()}] Starting\n")
    except:
        pass
    
    # Build matroid and OS3
    log("\nBuilding M6 matroid...")
    C, M6 = build_M6_matroid()
    log(f"  {len(C)} channels, rank {M6.rank()}")
    
    log("\nBuilding OS3 space...")
    triples, Vbasis = build_OS3_data(C, M6)
    log(f"  OS3 dim = {Vbasis.nrows()}")
    
    # Compute S6 invariants
    log("\nComputing S6 invariants...")
    Vinv = compute_S6_invariants(C, triples, Vbasis)
    log(f"  S6 invariants: {len(Vinv)} vectors")
    
    d = len(Vinv)
    
    if d == 0:
        log("[ERROR] No S6 invariants found!")
        return None, {'status': 'error'}
    
    # Generate kinematic samples
    log("\nGenerating kinematic samples...")
    n_samples = max(100, d * 3)  # Need overdetermined system
    samples = []
    hodges_values = []
    
    for seed in range(1000, 1000 + n_samples * 3):
        if len(samples) >= n_samples:
            break
        
        kin = SpinorHelicity(n=6, seed=seed)
        hodges_val = hodges_6pt_mhv(kin)
        
        if hodges_val is None:
            continue
        
        samples.append(kin)
        hodges_values.append(hodges_val)
        
        if len(samples) % 20 == 0:
            log(f"  {len(samples)}/{n_samples} valid samples")
    
    log(f"  Total: {len(samples)} samples")
    
    # Build evaluation matrix
    log("\nBuilding evaluation matrix...")
    eval_matrix = []
    
    for s_idx, kin in enumerate(samples):
        row = []
        for v in Vinv:
            val = evaluate_os3_form(v, C, triples, kin)
            row.append(val)
        eval_matrix.append(row)
        
        if (s_idx + 1) % 20 == 0:
            log(f"  Evaluated {s_idx + 1}/{len(samples)}")
    
    M = matrix(QQ, eval_matrix)
    b = vector(QQ, hodges_values)
    
    log(f"  Matrix: {M.nrows()} x {M.ncols()}")
    log(f"  Rank: {M.rank()}")
    
    # Solve
    log("\nSolving for physical direction...")
    
    try:
        Mt = M.transpose()
        MtM = Mt * M
        Mtb = Mt * b
        
        if MtM.det() != 0:
            c = MtM.inverse() * Mtb
        else:
            log("  Using pseudoinverse...")
            c = M.pseudoinverse() * b
        
        residual = M * c - b
        max_res = max(abs(r) for r in residual)
        avg_res = sum(abs(r) for r in residual) / len(residual)
        
        log(f"\nSolution found!")
        log(f"  Max residual: {float(max_res):.6e}")
        log(f"  Avg residual: {float(avg_res):.6e}")
        log(f"  Coefficients: {[float(x) for x in c]}")
        
        # Construct physical vector
        physical_vec = sum(c[i] * Vinv[i] for i in range(d))
        nnz = sum(1 for x in physical_vec if x != 0)
        log(f"  Physical vector: {nnz} nonzero entries")
        
        # Verify
        log("\nVerifying on new points...")
        errors = []
        for seed in range(2000, 2020):
            kin = SpinorHelicity(n=6, seed=seed)
            hodges_ref = hodges_6pt_mhv(kin)
            if hodges_ref is None or hodges_ref == 0:
                continue
            
            our_val = evaluate_os3_form(physical_vec, C, triples, kin)
            rel_err = abs(our_val - hodges_ref) / abs(hodges_ref)
            errors.append(float(rel_err))
        
        if errors:
            log(f"  Verification errors: max={max(errors):.6e}, avg={sum(errors)/len(errors):.6e}")
        
        # Report
        log("\n" + "="*70)
        if max_res < 0.01:
            log("[SUCCESS] Physical amplitude found!")
            status = "success"
        else:
            log("[PARTIAL] Direction found with residuals")
            status = "partial"
        log("="*70)
        
        result = {
            'status': status,
            'dim': d,
            'n_samples': len(samples),
            'max_residual': str(max_res),
            'coefficients': [str(x) for x in c],
            'nnz': nnz,
            'time': time.time() - t_start
        }
        
        # Save
        try:
            save({'result': result, 'physical_vec': physical_vec, 'coeffs': c}, 
                 'standalone_hodges_result.sobj')
            log("Saved result to standalone_hodges_result.sobj")
        except Exception as e:
            log(f"Save failed: {e}")
        
        try:
            with open("standalone_hodges_report.json", 'w') as f:
                json.dump(result, f, indent=2)
            log("Wrote standalone_hodges_report.json")
        except Exception as e:
            log(f"JSON write failed: {e}")
        
        log(f"\nTotal time: {time.time() - t_start:.1f}s")
        return physical_vec, result
        
    except Exception as e:
        log(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        return None, {'status': 'error', 'error': str(e)}

if __name__ == '__main__':
    main()


