#!/usr/bin/env sage
# =============================================================================
# DIRECT HODGES PROJECTION - Physics Breakthrough Strategy
# =============================================================================
#
# KEY INSIGHT FROM PREVIOUS RUNS:
# - S6 invariants (dim=2) become EMPTY under boundary constraints
# - This means gravity amplitude is NOT purely S6-invariant in OS3
# - BUT: S3xS3 after 1 boundary gives dim=48 (non-empty!)
#
# NEW STRATEGY:
# 1. Use S3xS3 invariants (dim=58) or after 1 boundary (dim=48)
# 2. Skip intersection - it makes things empty
# 3. Directly project onto Hodges amplitude using numerical evaluation
# 4. This finds the unique physical direction without needing intersection
#
# This is the "Phase C" approach from factorization_dim8_constraints_report.txt:
# "Match Hodges / CHY / KLT at points" projection to explicitly isolate the physical line
# =============================================================================

from sage.all import *
import numpy as np
import time
import os
import json
import gc
from itertools import combinations

# =============================================================================
# CONFIGURATION
# =============================================================================
DIAG = True
LOG_FILE = "hodges_projection.log"
CHECKPOINT_DIR = "hodges_checkpoints"

os.makedirs(CHECKPOINT_DIR, exist_ok=True)

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

def save_checkpoint(name, data):
    path = os.path.join(CHECKPOINT_DIR, f"{name}.sobj")
    try:
        save(data, path)
        log(f"  [CKPT] Saved: {name}")
    except Exception as e:
        log(f"  [CKPT] Failed: {name} - {e}")

def load_checkpoint(name):
    path = os.path.join(CHECKPOINT_DIR, f"{name}.sobj")
    if os.path.exists(path):
        try:
            return load(path)
        except:
            pass
    return None

# =============================================================================
# 4D SPINOR-HELICITY KINEMATICS
# =============================================================================

class SpinorHelicity:
    """
    4D spinor-helicity kinematics for n massless particles.
    
    Conventions:
    - p_i^{alpha dot{alpha}} = lambda_i^alpha * tilde_lambda_i^{dot{alpha}}
    - <ij> = epsilon_{alpha beta} lambda_i^alpha lambda_j^beta
    - [ij] = epsilon_{dot{alpha} dot{beta}} tilde_lambda_i^{dot{alpha}} tilde_lambda_j^{dot{beta}}
    - s_{ij} = <ij>[ji] = 2 p_i . p_j
    """
    
    def __init__(self, n=6, seed=42):
        self.n = n
        self.seed = seed
        np.random.seed(seed)
        
        # Generate random spinors with small integer components
        # Using rationals for exact arithmetic
        self.lambdas = []
        self.tilde_lambdas = []
        
        for i in range(n):
            # Random 2-component spinor
            lam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            tlam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            
            # Avoid zero spinors
            while lam[0] == 0 and lam[1] == 0:
                lam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            while tlam[0] == 0 and tlam[1] == 0:
                tlam = vector(QQ, [QQ(np.random.randint(-5, 6)), QQ(np.random.randint(-5, 6))])
            
            self.lambdas.append(lam)
            self.tilde_lambdas.append(tlam)
        
        # Precompute brackets and Mandelstams
        self._compute_brackets()
        self._compute_mandelstams()
    
    def _compute_brackets(self):
        """Compute all angle and square brackets."""
        n = self.n
        self.angle = {}  # <ij>
        self.square = {}  # [ij]
        
        for i in range(n):
            for j in range(n):
                # <ij> = lambda_i^1 * lambda_j^2 - lambda_i^2 * lambda_j^1
                self.angle[(i, j)] = self.lambdas[i][0] * self.lambdas[j][1] - self.lambdas[i][1] * self.lambdas[j][0]
                # [ij] = tilde_lambda_i^1 * tilde_lambda_j^2 - tilde_lambda_i^2 * tilde_lambda_j^1
                self.square[(i, j)] = self.tilde_lambdas[i][0] * self.tilde_lambdas[j][1] - self.tilde_lambdas[i][1] * self.tilde_lambdas[j][0]
    
    def _compute_mandelstams(self):
        """Compute all Mandelstam invariants."""
        n = self.n
        self.sij = {}  # s_{ij} for 2-particle
        self.sijk = {}  # s_{ijk} for 3-particle
        
        # 2-particle: s_{ij} = <ij>[ji]
        for i in range(n):
            for j in range(i+1, n):
                self.sij[(i+1, j+1)] = self.angle[(i, j)] * self.square[(j, i)]
        
        # 3-particle: s_{ijk} = (p_i + p_j + p_k)^2 = s_{ij} + s_{ik} + s_{jk}
        for triple in combinations(range(1, n+1), 3):
            i, j, k = triple
            val = self.sij.get((min(i,j), max(i,j)), QQ(0))
            val += self.sij.get((min(i,k), max(i,k)), QQ(0))
            val += self.sij.get((min(j,k), max(j,k)), QQ(0))
            self.sijk[triple] = val
    
    def get_sij(self, i, j):
        """Get s_{ij} (1-indexed)."""
        if i > j:
            i, j = j, i
        return self.sij.get((i, j), QQ(0))
    
    def get_sijk(self, i, j, k):
        """Get s_{ijk} (1-indexed)."""
        triple = tuple(sorted([i, j, k]))
        return self.sijk.get(triple, QQ(0))

# =============================================================================
# HODGES FORMULA FOR 6-POINT MHV GRAVITY
# =============================================================================

def hodges_6pt_mhv(kin):
    """
    Compute Hodges formula for 6-point MHV gravity amplitude.
    
    For MHV with particles 1,2 negative helicity, rest positive:
    M_6 = det'(Phi) / (<12><23><34><45><56><61>)
    
    where Phi is the reduced (n-2)x(n-2) matrix (delete rows/cols for 1 and n).
    
    Phi_{ij} = [ij]/<ij> for i != j (i,j in {2,3,4,5})
    Phi_{ii} = -sum_{k != i,1,6} [ik]<1k><6k>/(<ik><1i><6i>)
    """
    n = kin.n
    
    # Build Phi matrix for indices 2,3,4,5 (0-indexed: 1,2,3,4)
    indices = [1, 2, 3, 4]  # 0-indexed
    d = len(indices)
    
    Phi = matrix(QQ, d, d)
    
    for ii, i in enumerate(indices):
        for jj, j in enumerate(indices):
            if ii == jj:
                # Diagonal element
                diag_sum = QQ(0)
                for k in range(n):
                    if k in [i, 0, 5]:  # Skip i, particle 1 (0), particle 6 (5)
                        continue
                    
                    ik_angle = kin.angle[(i, k)]
                    i1_angle = kin.angle[(i, 0)]  # <i,1>
                    i6_angle = kin.angle[(i, 5)]  # <i,6>
                    
                    if ik_angle == 0 or i1_angle == 0 or i6_angle == 0:
                        continue
                    
                    # [ik]<1k><6k>/(<ik><1i><6i>)
                    contrib = kin.square[(i, k)] * kin.angle[(0, k)] * kin.angle[(5, k)]
                    contrib = contrib / (ik_angle * i1_angle * i6_angle)
                    diag_sum -= contrib
                
                Phi[ii, jj] = diag_sum
            else:
                # Off-diagonal: [ij]/<ij>
                ij_angle = kin.angle[(i, j)]
                if ij_angle == 0:
                    return None  # Singular
                Phi[ii, jj] = kin.square[(i, j)] / ij_angle
    
    # Compute det(Phi)
    try:
        det_Phi = Phi.det()
    except:
        return None
    
    # Denominator: <12><23><34><45><56><61>
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

# =============================================================================
# EVALUATE OS3 FORM ON KINEMATICS
# =============================================================================

def evaluate_os3_form(coeffs, C, triples, kin):
    """
    Evaluate an OS3 differential form on given kinematics.
    
    The form is: sum_{(i,j,k) in triples} c_{ijk} * d(log s_i) ^ d(log s_j) ^ d(log s_k)
    
    At a kinematic point, this becomes a rational function of Mandelstams.
    We evaluate it as: sum c_{ijk} / (s_i * s_j * s_k)
    
    This is a simplification - the full evaluation requires more careful treatment
    of the differential form structure.
    """
    result = QQ(0)
    
    for t_idx, (i, j, k) in enumerate(triples):
        if t_idx >= len(coeffs):
            break
        
        c = coeffs[t_idx]
        if c == 0:
            continue
        
        # Get channel values
        ch_i = C[i]
        ch_j = C[j]
        ch_k = C[k]
        
        # Evaluate each channel
        def get_channel_value(ch):
            if len(ch) == 2:
                return kin.get_sij(ch[0], ch[1])
            else:
                return kin.get_sijk(ch[0], ch[1], ch[2])
        
        si = get_channel_value(ch_i)
        sj = get_channel_value(ch_j)
        sk = get_channel_value(ch_k)
        
        if si == 0 or sj == 0 or sk == 0:
            continue
        
        result += c / (si * sj * sk)
    
    return result

# =============================================================================
# MAIN: DIRECT HODGES PROJECTION
# =============================================================================

def main():
    """Direct projection to Hodges amplitude."""
    
    log("\n" + "="*70)
    log("DIRECT HODGES PROJECTION - Physics Breakthrough")
    log("="*70)
    log("Strategy: Skip boundary intersection (causes empty space)")
    log("Instead: Project candidate space directly onto Hodges amplitude")
    log("="*70)
    
    t_start = time.time()
    
    # Clear log
    try:
        with open(LOG_FILE, 'w') as f:
            f.write(f"[{ts()}] Starting direct Hodges projection\n")
    except:
        pass
    
    # Phase 1: Load infrastructure from 54.sage
    log("\n" + "="*70)
    log("PHASE 1: LOAD MATROID AND OS3 INFRASTRUCTURE")
    log("="*70)
    
    # Load 54.sage but override settings
    log("Loading 54.sage infrastructure...")
    
    # Set overrides
    global INTERSECTION_MODE, MULTI_STRATEGY_SEARCH
    INTERSECTION_MODE = False
    MULTI_STRATEGY_SEARCH = False
    
    load('54.sage')
    
    log("Infrastructure loaded.")
    
    # Build matroid and OS3
    log("\nBuilding M6 matroid...")
    C, M6 = build_M6_matroid()
    log(f"  {len(C)} channels, rank {M6.rank()}")
    
    log("\nBuilding OS3 space...")
    triples, Vbasis = build_OS3_data(C, M6)
    log(f"  OS3 dim = {Vbasis.nrows()}")
    
    # Phase 2: Compute S3xS3 invariants (larger than S6)
    log("\n" + "="*70)
    log("PHASE 2: COMPUTE S3xS3 INVARIANTS")
    log("="*70)
    
    S3xS3_gens = [
        SymmetricGroup(6)((1, 2)),
        SymmetricGroup(6)((2, 3)),
        SymmetricGroup(6)((4, 5)),
        SymmetricGroup(6)((5, 6)),
    ]
    
    Vinv = compute_invariants_from_generators(
        C, triples, Vbasis, S3xS3_gens, "invariants_S3xS3.sobj"
    )
    log(f"  S3xS3 invariants: {len(Vinv)} vectors")
    
    # Phase 3: Generate kinematic sample points
    log("\n" + "="*70)
    log("PHASE 3: GENERATE KINEMATIC SAMPLES")
    log("="*70)
    
    n_samples = 100  # Use many samples for robust projection
    samples = []
    hodges_values = []
    
    valid_count = 0
    for seed in range(1000, 1000 + n_samples * 2):  # Try more seeds to get enough valid ones
        if valid_count >= n_samples:
            break
        
        kin = SpinorHelicity(n=6, seed=seed)
        hodges_val = hodges_6pt_mhv(kin)
        
        if hodges_val is None:
            continue  # Skip singular kinematics
        
        samples.append(kin)
        hodges_values.append(hodges_val)
        valid_count += 1
        
        if valid_count % 20 == 0:
            log(f"  Generated {valid_count}/{n_samples} valid samples")
    
    log(f"  Total valid samples: {len(samples)}")
    
    if len(samples) < len(Vinv):
        log(f"[WARNING] Need at least {len(Vinv)} samples, only have {len(samples)}")
        log("  Trying to continue anyway...")
    
    # Phase 4: Build evaluation matrix
    log("\n" + "="*70)
    log("PHASE 4: BUILD EVALUATION MATRIX")
    log("="*70)
    
    d = len(Vinv)
    n_pts = len(samples)
    
    log(f"Building {n_pts} x {d} evaluation matrix...")
    
    eval_matrix = []
    
    for s_idx, kin in enumerate(samples):
        row = []
        for v_idx, v in enumerate(Vinv):
            # Evaluate basis vector v at kinematic point kin
            val = evaluate_os3_form(v, C, triples, kin)
            row.append(val)
        eval_matrix.append(row)
        
        if (s_idx + 1) % 20 == 0:
            log(f"  Evaluated {s_idx + 1}/{n_pts} samples")
    
    M = matrix(QQ, eval_matrix)
    b = vector(QQ, hodges_values)
    
    log(f"  Matrix size: {M.nrows()} x {M.ncols()}")
    log(f"  Matrix rank: {M.rank()}")
    
    # Phase 5: Solve for physical direction
    log("\n" + "="*70)
    log("PHASE 5: SOLVE FOR PHYSICAL DIRECTION")
    log("="*70)
    
    try:
        # Try to solve M * c = b
        log("Attempting to solve M * c = b...")
        
        # Check if system is consistent
        augmented = M.augment(b.column())
        rank_M = M.rank()
        rank_aug = augmented.rank()
        
        log(f"  rank(M) = {rank_M}, rank([M|b]) = {rank_aug}")
        
        if rank_M != rank_aug:
            log("[WARNING] System may be inconsistent")
        
        # Use least squares if overdetermined
        if M.nrows() > M.ncols():
            log("Using least squares (overdetermined system)...")
            Mt = M.transpose()
            MtM = Mt * M
            Mtb = Mt * b
            
            if MtM.det() != 0:
                c = MtM.inverse() * Mtb
                log("  Solved via normal equations")
            else:
                log("  MtM is singular, using pseudoinverse...")
                # SVD-based pseudoinverse
                c = M.pseudoinverse() * b
        else:
            c = M.solve_right(b)
            log("  Solved exactly")
        
        # Check residual
        residual = M * c - b
        max_res = max(abs(r) for r in residual)
        avg_res = sum(abs(r) for r in residual) / len(residual)
        
        log(f"\nSolution found!")
        log(f"  Max residual: {float(max_res):.6e}")
        log(f"  Avg residual: {float(avg_res):.6e}")
        
        # Construct physical vector in OS3 space
        physical_vec = sum(c[i] * Vinv[i] for i in range(d))
        
        # Count nonzeros
        nnz = sum(1 for x in physical_vec if x != 0)
        log(f"  Physical vector: {nnz} nonzero entries")
        
        # Verify on a few more random points
        log("\nVerifying on additional random points...")
        n_verify = 10
        verify_errors = []
        
        for seed in range(2000, 2000 + n_verify * 2):
            kin = SpinorHelicity(n=6, seed=seed)
            hodges_ref = hodges_6pt_mhv(kin)
            if hodges_ref is None:
                continue
            
            our_val = evaluate_os3_form(physical_vec, C, triples, kin)
            if hodges_ref != 0:
                rel_error = abs(our_val - hodges_ref) / abs(hodges_ref)
                verify_errors.append(float(rel_error))
            
            if len(verify_errors) >= n_verify:
                break
        
        if verify_errors:
            max_verify_error = max(verify_errors)
            avg_verify_error = sum(verify_errors) / len(verify_errors)
            log(f"  Verification max error: {max_verify_error:.6e}")
            log(f"  Verification avg error: {avg_verify_error:.6e}")
        
        # Save results
        log("\n" + "="*70)
        log("RESULTS")
        log("="*70)
        
        if max_res < QQ(1)/QQ(100):
            log("[SUCCESS] Physical amplitude direction found!")
            status = "success"
        else:
            log("[PARTIAL] Direction found but with significant residuals")
            status = "partial"
        
        result = {
            'status': status,
            'candidate_dim': d,
            'n_samples': len(samples),
            'max_residual': str(max_res),
            'avg_residual': str(avg_res),
            'nnz': nnz,
            'coefficients': [str(x) for x in c],
            'verification_errors': verify_errors,
            'timestamp': time.strftime("%Y-%m-%d %H:%M:%S"),
            'total_time': time.time() - t_start
        }
        
        save_checkpoint("hodges_projection_result", result)
        save_checkpoint("physical_vector", physical_vec)
        save_checkpoint("coefficients", c)
        
        # Write JSON report
        try:
            with open("hodges_projection_report.json", 'w') as f:
                json.dump(result, f, indent=2)
            log("Wrote hodges_projection_report.json")
        except Exception as e:
            log(f"Failed to write JSON: {e}")
        
        log(f"\nTotal time: {time.time() - t_start:.1f}s")
        
        return physical_vec, result
        
    except Exception as e:
        log(f"[ERROR] Projection failed: {e}")
        import traceback
        traceback.print_exc()
        return None, {'status': 'error', 'error': str(e)}

if __name__ == '__main__':
    main()










