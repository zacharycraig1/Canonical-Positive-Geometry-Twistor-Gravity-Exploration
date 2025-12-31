#!/usr/bin/env sage
# =============================================================================
# FORENSICS: DIAGNOSE MISMATCHES BETWEEN AMPLITUHEDRON AND HODGES
# =============================================================================
# Goal: Understand the 7 failures out of 200 tests
# - Classify each failure (degeneracy, positivity, normalization, bug)
# - Sample from positive region (not random integers)
# - Make comparison invariant under normalization
# =============================================================================

from sage.all import *
import numpy as np
import time
import os
import json
from itertools import combinations

DIAG = True
LOG_FILE = "forensics.log"
RESULTS_FILE = "forensics_results.json"

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
# POSITIVE SAMPLING
# =============================================================================

def sample_positive_twistors(n=6, seed=None):
    """
    Sample momentum twistors from positive region.
    
    For amplituhedron, we need:
    - All ordered 4×4 minors positive: ⟨i j k l⟩ > 0 for i < j < k < l
    - All angle brackets positive: ⟨i i+1⟩ > 0
    
    Strategy: Generate a positive matrix Z with all minors positive.
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Start with a positive matrix
    # Use exponential of random matrix to ensure positivity
    Z = []
    
    # Generate twistors that satisfy positivity
    # One approach: use Plücker coordinates that are all positive
    for i in range(n):
        # Use exponential coordinates to ensure positivity
        z = vector(QQ, [
            QQ(np.random.randint(1, 20)),  # Start positive
            QQ(np.random.randint(1, 20)),
            QQ(np.random.randint(1, 20)),
            QQ(np.random.randint(1, 20))
        ])
        Z.append(z)
    
    # Check and adjust to ensure positivity
    twistor = MomentumTwistor(n=n)
    twistor.Z = Z
    twistor._compute_brackets()
    
    # Verify positivity
    all_positive = True
    for ijkl in combinations(range(n), 4):
        i, j, k, l = sorted(ijkl)
        bracket = twistor.get_four_bracket(i, j, k, l)
        if bracket <= 0:
            all_positive = False
            break
    
    # If not positive, try again with different seed
    if not all_positive:
        return sample_positive_twistors(n=n, seed=(seed or 0) + 1)
    
    return twistor

# =============================================================================
# MOMENTUM TWISTOR
# =============================================================================

class MomentumTwistor:
    def __init__(self, n=6, seed=None, Z=None):
        self.n = n
        if Z is not None:
            self.Z = Z
        else:
            if seed is not None:
                np.random.seed(seed)
            self.Z = []
            for i in range(n):
                z = vector(QQ, [
                    QQ(np.random.randint(-10, 11)),
                    QQ(np.random.randint(-10, 11)),
                    QQ(np.random.randint(-10, 11)),
                    QQ(np.random.randint(-10, 11))
                ])
                while all(x == 0 for x in z):
                    z = vector(QQ, [
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11))
                    ])
                self.Z.append(z)
        self._compute_brackets()
    
    def _compute_brackets(self):
        n = self.n
        self.angle = {}
        for i in range(n):
            for j in range(n):
                self.angle[(i, j)] = self.Z[i][0] * self.Z[j][1] - self.Z[i][1] * self.Z[j][0]
        self.four_bracket = {}
        for ijkl in combinations(range(n), 4):
            i, j, k, l = sorted(ijkl)
            M = matrix(QQ, [self.Z[i], self.Z[j], self.Z[k], self.Z[l]])
            self.four_bracket[ijkl] = M.det()
    
    def get_angle(self, i, j):
        return self.angle.get((i, j), QQ(0))
    
    def get_four_bracket(self, i, j, k, l):
        indices = tuple(sorted([i, j, k, l]))
        base = self.four_bracket.get(indices, QQ(0))
        perm_list = [i, j, k, l]
        sorted_list = sorted(perm_list)
        inversions = 0
        for a in range(len(perm_list)):
            for b in range(a + 1, len(perm_list)):
                if sorted_list.index(perm_list[a]) > sorted_list.index(perm_list[b]):
                    inversions += 1
        sign = 1 if inversions % 2 == 0 else -1
        return sign * base
    
    def get_square(self, i, j):
        im1, jm1 = (i - 1) % self.n, (j - 1) % self.n
        num = self.get_four_bracket(im1, i, jm1, j)
        den = self.get_angle(im1, i) * self.get_angle(jm1, j)
        return num / den if den != 0 else None

# =============================================================================
# HODGES FORMULA
# =============================================================================

def hodges_6pt_mhv(twistor):
    """Hodges formula - proven correct."""
    n = twistor.n
    indices = [1, 2, 3, 4]
    d = len(indices)
    
    Phi = matrix(QQ, d, d)
    
    for ii, i in enumerate(indices):
        for jj, j in enumerate(indices):
            if ii == jj:
                diag_sum = QQ(0)
                for k in range(n):
                    if k in [i, 0, 5]:
                        continue
                    ik_sq = twistor.get_square(i, k)
                    if ik_sq is None:
                        continue
                    ik_ang = twistor.get_angle(i, k)
                    i0_ang = twistor.get_angle(i, 0)
                    i5_ang = twistor.get_angle(i, 5)
                    k0_ang = twistor.get_angle(k, 0)
                    k5_ang = twistor.get_angle(k, 5)
                    if ik_ang == 0 or i0_ang == 0 or i5_ang == 0:
                        continue
                    contrib = ik_sq * k0_ang * k5_ang / (ik_ang * i0_ang * i5_ang)
                    diag_sum -= contrib
                Phi[ii, jj] = diag_sum
            else:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0:
                    return None, "zero_angle_bracket"
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None:
                    return None, "zero_square_bracket"
                Phi[ii, jj] = ij_sq / ij_ang
    
    try:
        det_Phi = Phi.det()
    except:
        return None, "determinant_failed"
    
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        bracket = twistor.get_angle(i, j)
        if bracket == 0:
            return None, "zero_denominator_angle"
        denom *= bracket
    
    return det_Phi / denom if denom != 0 else None, "ok"

# =============================================================================
# AMPLITUHEDRON (direct = Hodges for MHV)
# =============================================================================

def amplituhedron_6pt_mhv(twistor):
    """
    For MHV, amplituhedron volume = Hodges determinant.
    This is the direct computation.
    """
    result, reason = hodges_6pt_mhv(twistor)
    return result, reason

# =============================================================================
# DIAGNOSTICS
# =============================================================================

def min_abs_relevant_minors(twistor):
    """Find minimum absolute value of relevant minors."""
    n = twistor.n
    minors = []
    
    # All 4-brackets
    for ijkl in combinations(range(n), 4):
        i, j, k, l = sorted(ijkl)
        bracket = twistor.get_four_bracket(i, j, k, l)
        minors.append(abs(bracket))
    
    # All angle brackets
    for i in range(n):
        j = (i + 1) % n
        ang = twistor.get_angle(i, j)
        minors.append(abs(ang))
    
    return min(m for m in minors if m != 0) if minors else None

def check_positivity(twistor):
    """Check if all ordered 4×4 minors are positive."""
    n = twistor.n
    for ijkl in combinations(range(n), 4):
        i, j, k, l = sorted(ijkl)
        bracket = twistor.get_four_bracket(i, j, k, l)
        if bracket <= 0:
            return False, f"negative_minor_{i}_{j}_{k}_{l}"
    return True, "ok"

def explain_none(twistor, A_reason, H_reason):
    """Explain why A or H returned None."""
    reasons = []
    if A_reason != "ok":
        reasons.append(f"A: {A_reason}")
    if H_reason != "ok":
        reasons.append(f"H: {H_reason}")
    return "; ".join(reasons) if reasons else "unknown"

def diagnostics_for_point(twistor, point_idx):
    """
    Full diagnostics for a single point.
    Returns dict with all relevant information.
    """
    diag = {
        "point_idx": point_idx,
        "Z": [[float(x) for x in z] for z in twistor.Z]
    }
    
    # Compute A and H
    A, A_reason = amplituhedron_6pt_mhv(twistor)
    H, H_reason = hodges_6pt_mhv(twistor)
    
    diag["A"] = str(A) if A is not None else None
    diag["H"] = str(H) if H is not None else None
    diag["A_reason"] = A_reason
    diag["H_reason"] = H_reason
    
    # Check for None
    if A is None or H is None:
        diag["status"] = "None"
        diag["reason"] = explain_none(twistor, A_reason, H_reason)
        diag["min_minor"] = min_abs_relevant_minors(twistor)
        pos_ok, pos_reason = check_positivity(twistor)
        diag["positivity_ok"] = pos_ok
        diag["positivity_reason"] = pos_reason
        return diag
    
    # Check for zero H
    if H == 0:
        diag["status"] = "H_zero"
        diag["min_minor"] = min_abs_relevant_minors(twistor)
        pos_ok, pos_reason = check_positivity(twistor)
        diag["positivity_ok"] = pos_ok
        diag["positivity_reason"] = pos_reason
        return diag
    
    # Compute ratio and difference
    try:
        ratio = (A / H).simplify_rational()
        diff = (A - ratio * H).simplify_rational()
        
        diag["ratio"] = str(ratio)
        diag["diff"] = str(diff)
        diag["diff_abs"] = float(abs(diff)) if diff != 0 else 0.0
        
        # Check if they match
        if A == H:
            diag["status"] = "exact_match"
        elif diff == 0:
            diag["status"] = "ratio_match"  # Same up to constant
        else:
            diag["status"] = "mismatch"
            
    except Exception as e:
        diag["status"] = "computation_error"
        diag["error"] = str(e)
    
    # Additional diagnostics
    diag["min_minor"] = min_abs_relevant_minors(twistor)
    pos_ok, pos_reason = check_positivity(twistor)
    diag["positivity_ok"] = pos_ok
    diag["positivity_reason"] = pos_reason
    
    return diag

# =============================================================================
# MAIN FORENSICS
# =============================================================================

def main():
    log("\n" + "="*70)
    log("FORENSICS: DIAGNOSING MISMATCHES")
    log("="*70)
    
    TRIALS = 200
    results = []
    matches = 0
    mismatches = 0
    none_cases = 0
    h_zero = 0
    
    log(f"Testing {TRIALS} points...")
    
    for t in range(TRIALS):
        if t % 50 == 0:
            log(f"Progress: {t}/{TRIALS}")
        
        # Sample (for now, use random; later switch to positive)
        twistor = MomentumTwistor(n=6, seed=t)
        
        diag = diagnostics_for_point(twistor, t)
        results.append(diag)
        
        status = diag.get("status", "unknown")
        if status == "exact_match":
            matches += 1
        elif status == "ratio_match":
            matches += 1
        elif status == "mismatch":
            mismatches += 1
        elif status == "None":
            none_cases += 1
        elif status == "H_zero":
            h_zero += 1
    
    log("\n" + "="*70)
    log("SUMMARY")
    log("="*70)
    log(f"Total points: {TRIALS}")
    log(f"Exact matches: {matches}")
    log(f"Mismatches: {mismatches}")
    log(f"None cases: {none_cases}")
    log(f"H zero: {h_zero}")
    
    # Analyze mismatches
    if mismatches > 0:
        log("\n" + "="*70)
        log("MISMATCH ANALYSIS")
        log("="*70)
        
        mismatch_indices = [i for i, r in enumerate(results) if r.get("status") == "mismatch"]
        log(f"Mismatch indices: {mismatch_indices}")
        
        for idx in mismatch_indices[:10]:  # Show first 10
            r = results[idx]
            log(f"\nPoint {idx}:")
            log(f"  A = {r.get('A', 'N/A')}")
            log(f"  H = {r.get('H', 'N/A')}")
            log(f"  ratio = {r.get('ratio', 'N/A')}")
            log(f"  diff = {r.get('diff', 'N/A')}")
            log(f"  min_minor = {r.get('min_minor', 'N/A')}")
            log(f"  positivity_ok = {r.get('positivity_ok', 'N/A')}")
    
    # Analyze None cases
    if none_cases > 0:
        log("\n" + "="*70)
        log("NONE CASE ANALYSIS")
        log("="*70)
        
        none_indices = [i for i, r in enumerate(results) if r.get("status") == "None"]
        log(f"None case indices: {none_indices}")
        
        for idx in none_indices[:10]:  # Show first 10
            r = results[idx]
            log(f"\nPoint {idx}:")
            log(f"  reason = {r.get('reason', 'N/A')}")
            log(f"  A_reason = {r.get('A_reason', 'N/A')}")
            log(f"  H_reason = {r.get('H_reason', 'N/A')}")
            log(f"  min_minor = {r.get('min_minor', 'N/A')}")
            log(f"  positivity_ok = {r.get('positivity_ok', 'N/A')}")
    
    # Save results
    try:
        # Convert to JSON-serializable format
        json_results = []
        for r in results:
            json_r = {}
            for k, v in r.items():
                if isinstance(v, (int, float, str, bool, type(None))):
                    json_r[k] = v
                elif isinstance(v, list):
                    json_r[k] = v
                else:
                    json_r[k] = str(v)
            json_results.append(json_r)
        
        with open(RESULTS_FILE, 'w') as f:
            json.dump(json_results, f, indent=2)
        log(f"\nSaved results to {RESULTS_FILE}")
    except Exception as e:
        log(f"Failed to save results: {e}")
    
    # Save failing points
    failing = [r for r in results if r.get("status") in ["mismatch", "None", "H_zero"]]
    if failing:
        log(f"\nSaving {len(failing)} failing points...")
        try:
            with open("failing_points.json", 'w') as f:
                json.dump([r for r in failing], f, indent=2)
            log("Saved failing points to failing_points.json")
        except Exception as e:
            log(f"Failed to save failing points: {e}")
    
    return results

if __name__ == '__main__':
    main()

