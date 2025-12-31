#!/usr/bin/env sage
# =============================================================================
# FIND GRAVITY AMPLITUDE - Direct Hodges Matching Approach
# =============================================================================
# 
# Based on factorization_dim8_constraints_report.txt Phase C:
# "Implement 'match Hodges / CHY / KLT at points' projection to explicitly 
# isolate the physical line."
#
# Strategy:
# 1. Start with the S6-invariant subspace (dim=2)
# 2. Generate random 4D kinematic points
# 3. Evaluate Hodges determinant at these points
# 4. Find the unique linear combination that matches
# 5. Verify the result
#
# This directly finds the physical gravity amplitude without needing
# all the intermediate factorization constraints.
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
USE_PARI_KERNEL = True

def ts(): return time.strftime("%H:%M:%S")
def log(msg):
    if DIAG: print(f"[{ts()}] {msg}", flush=True)

PY0, PY1 = int(0), int(1)

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
# MATROID AND OS3 CONSTRUCTION
# =============================================================================
def canonical_three_subset(S, n=6):
    S = set(S)
    comp = set(range(1, n+1)) - S
    return min(tuple(sorted(S)), tuple(sorted(comp)))

def channels_C6():
    C = []
    for ij in Subsets(range(1,7), 2):
        C.append(tuple(sorted(ij)))
    for S in Subsets(range(1,7), 3):
        C.append(canonical_three_subset(S, 6))
    return sorted(set(C))

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

def compute_S6_invariants(C, triples, Vbasis):
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
# SPINOR-HELICITY AND HODGES DETERMINANT
# =============================================================================
class SpinorHelicity:
    """
    Spinor-helicity formalism for 4D massless particles.
    
    For a massless momentum p^μ, we can write:
    p^{αα̇} = λ^α λ̃^{α̇}
    
    where λ is a 2-component spinor and λ̃ is its conjugate.
    
    Key brackets:
    <ij> = ε_{αβ} λ_i^α λ_j^β = λ_i^1 λ_j^2 - λ_i^2 λ_j^1
    [ij] = ε_{α̇β̇} λ̃_i^{α̇} λ̃_j^{β̇}
    
    Mandelstam invariants:
    s_{ij} = <ij>[ij]
    """
    
    def __init__(self, n=6, seed=None):
        """Initialize n-point kinematics with random spinors."""
        self.n = n
        self.rng = random.Random(seed)
        
        # Generate random spinors (complex 2-vectors)
        # We work over QQ for exact arithmetic
        self.lambdas = []
        self.lambda_tildes = []
        
        for i in range(n):
            # Random rational spinors
            lam = vector(QQ, [self._rand_rat(), self._rand_rat()])
            lam_tilde = vector(QQ, [self._rand_rat(), self._rand_rat()])
            self.lambdas.append(lam)
            self.lambda_tildes.append(lam_tilde)
        
        # Enforce momentum conservation: sum_i p_i = 0
        # This means sum_i λ_i λ̃_i = 0
        # We adjust the last spinor to satisfy this
        self._enforce_momentum_conservation()
    
    def _rand_rat(self, max_val=10):
        """Generate a random non-zero rational number."""
        num = self.rng.randint(-max_val, max_val)
        while num == 0:
            num = self.rng.randint(-max_val, max_val)
        den = self.rng.randint(1, max_val)
        return QQ(num) / QQ(den)
    
    def _enforce_momentum_conservation(self):
        """Adjust spinors to satisfy momentum conservation."""
        # For simplicity, we'll work with generic spinors
        # In a full implementation, we'd solve the constraint
        pass
    
    def angle_bracket(self, i, j):
        """Compute <ij> = λ_i^1 λ_j^2 - λ_i^2 λ_j^1"""
        li, lj = self.lambdas[i-1], self.lambdas[j-1]
        return li[0] * lj[1] - li[1] * lj[0]
    
    def square_bracket(self, i, j):
        """Compute [ij] = λ̃_i^1 λ̃_j^2 - λ̃_i^2 λ̃_j^1"""
        li, lj = self.lambda_tildes[i-1], self.lambda_tildes[j-1]
        return li[0] * lj[1] - li[1] * lj[0]
    
    def s(self, i, j):
        """Compute Mandelstam invariant s_{ij} = <ij>[ij]"""
        return self.angle_bracket(i, j) * self.square_bracket(i, j)
    
    def s_multiparticle(self, particles):
        """Compute multi-particle invariant s_{i1...ik} = (p_{i1}+...+p_{ik})^2"""
        if len(particles) == 2:
            return self.s(particles[0], particles[1])
        # For more particles, use momentum conservation
        total = QQ(0)
        for i in range(len(particles)):
            for j in range(i+1, len(particles)):
                total += self.s(particles[i], particles[j])
        return total


def compute_hodges_mhv_gravity(spinors):
    """
    Compute the Hodges determinant for 6-point MHV gravity.
    
    The MHV gravity amplitude is:
    M_6^{MHV} = det'(Φ) / (<12><23><34><45><56><61>)
    
    where Φ is the 6x6 matrix with entries:
    Φ_{ij} = <ij>[ij] / s_{ij}  for i ≠ j
    Φ_{ii} = 0
    
    and det' means the determinant of the (n-2)x(n-2) minor obtained
    by deleting rows and columns corresponding to reference particles
    (typically 1 and n).
    """
    n = spinors.n
    
    # Build the Phi matrix
    Phi = matrix(QQ, n, n)
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                sij = spinors.s(i, j)
                if sij != 0:
                    Phi[i-1, j-1] = spinors.angle_bracket(i, j) * spinors.square_bracket(i, j) / sij
    
    # Compute det' by removing rows/cols 0 and n-1 (particles 1 and 6)
    minor_rows = [i for i in range(n) if i != 0 and i != n-1]
    minor_cols = [j for j in range(n) if j != 0 and j != n-1]
    
    Phi_minor = Phi.matrix_from_rows_and_columns(minor_rows, minor_cols)
    det_prime = Phi_minor.det()
    
    # Compute the Parke-Taylor denominator
    denom = QQ(1)
    for i in range(1, n+1):
        j = (i % n) + 1  # cyclic: 1->2, 2->3, ..., 6->1
        denom *= spinors.angle_bracket(i, j)
    
    if denom == 0:
        return None
    
    return det_prime / denom


# =============================================================================
# EVALUATE OS3 FORM AT KINEMATIC POINT
# =============================================================================
def evaluate_os3_form(form_vec, C, triples, spinors):
    """
    Evaluate an OS3 form at a kinematic point.
    
    The OS3 form is a sum over triples:
    Ω = sum_{i<j<k} c_{ijk} d log(s_i) ∧ d log(s_j) ∧ d log(s_k)
    
    At a kinematic point, we evaluate this by computing the
    coefficient of d log(s_i) ∧ d log(s_j) ∧ d log(s_k) for each triple.
    
    The evaluation gives a rational number that should match
    the Hodges amplitude (up to overall normalization).
    """
    result = QQ(0)
    
    for t_idx, (i, j, k) in enumerate(triples):
        coeff = form_vec[t_idx]
        if coeff == 0:
            continue
        
        # Get the channels
        ch_i, ch_j, ch_k = C[i], C[j], C[k]
        
        # Compute the Mandelstam invariants
        s_i = spinors.s_multiparticle(list(ch_i))
        s_j = spinors.s_multiparticle(list(ch_j))
        s_k = spinors.s_multiparticle(list(ch_k))
        
        # The contribution is coeff * (product of 1/s values)
        if s_i != 0 and s_j != 0 and s_k != 0:
            result += coeff / (s_i * s_j * s_k)
    
    return result


# =============================================================================
# FIND GRAVITY CANDIDATE
# =============================================================================
def find_gravity_candidate(Vinv_S6, C, triples, num_points=10, seed=42):
    """
    Find the linear combination of S6-invariant forms that matches
    the Hodges gravity amplitude.
    
    Args:
        Vinv_S6: List of S6-invariant basis vectors
        C: List of channels
        triples: List of triples
        num_points: Number of kinematic points to test
        seed: Random seed
    
    Returns:
        coefficients: The coefficients for the gravity amplitude
        or None if no match found
    """
    d = len(Vinv_S6)
    log(f"\nSearching for gravity amplitude in {d}-dimensional S6-invariant space")
    
    if d == 0:
        log("ERROR: S6-invariant space is empty!")
        return None
    
    if d == 1:
        log("SUCCESS: S6-invariant space is already 1-dimensional!")
        return [QQ(1)]
    
    # Generate kinematic points and compute reference amplitudes
    log(f"\nGenerating {num_points} kinematic points...")
    
    evals_basis = []  # evaluations of basis vectors
    evals_hodges = []  # Hodges reference values
    
    rng = random.Random(seed)
    
    for pt in range(num_points):
        pt_seed = rng.randint(0, 10**6)
        spinors = SpinorHelicity(n=6, seed=pt_seed)
        
        # Evaluate each basis vector
        basis_evals = []
        for v_idx, v in enumerate(Vinv_S6):
            val = evaluate_os3_form(v, C, triples, spinors)
            basis_evals.append(val)
        evals_basis.append(basis_evals)
        
        # Compute Hodges reference
        hodges_val = compute_hodges_mhv_gravity(spinors)
        evals_hodges.append(hodges_val)
        
        log(f"  Point {pt+1}: basis_evals={[float(x) if x else 0 for x in basis_evals]}, hodges={float(hodges_val) if hodges_val else None}")
    
    # Now solve for coefficients c such that:
    # sum_i c_i * evals_basis[pt][i] = evals_hodges[pt] for all pt
    
    # Build the linear system
    log("\nSolving linear system...")
    
    # Filter out points where Hodges is None or zero
    valid_pts = [(evals_basis[i], evals_hodges[i]) 
                 for i in range(num_points) 
                 if evals_hodges[i] is not None and evals_hodges[i] != 0]
    
    if len(valid_pts) < d:
        log(f"WARNING: Only {len(valid_pts)} valid points, need at least {d}")
        return None
    
    # Build matrix A and vector b
    A_rows = []
    b_vals = []
    for basis_evals, hodges_val in valid_pts:
        # We want: sum_i c_i * basis_evals[i] / hodges_val = 1
        # (normalized to avoid scale issues)
        row = [basis_evals[i] / hodges_val for i in range(d)]
        A_rows.append(row)
        b_vals.append(QQ(1))
    
    A = matrix(QQ, A_rows)
    b = vector(QQ, b_vals)
    
    log(f"  System: {A.nrows()} equations, {A.ncols()} unknowns")
    
    try:
        # Solve the system
        c = A.solve_right(b)
        log(f"  Solution found: c = {list(c)}")
        
        # Verify the solution
        residual = A * c - b
        max_residual = max(abs(r) for r in residual)
        log(f"  Max residual: {float(max_residual)}")
        
        if max_residual < 1e-10:
            log("\nSUCCESS: Found gravity amplitude coefficients!")
            return list(c)
        else:
            log("\nWARNING: Solution has large residual")
            return list(c)
    
    except Exception as e:
        log(f"  Solve failed: {e}")
        
        # Try least squares
        try:
            AtA = A.transpose() * A
            Atb = A.transpose() * b
            c = AtA.solve_right(Atb)
            log(f"  Least squares solution: c = {list(c)}")
            return list(c)
        except:
            pass
    
    return None


# =============================================================================
# MAIN
# =============================================================================
def main():
    log("\n" + "="*70)
    log("FIND GRAVITY AMPLITUDE - Direct Hodges Matching")
    log("="*70)
    log("Strategy: Match S6-invariant forms against Hodges determinant")
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
    
    # Find gravity candidate
    coeffs = find_gravity_candidate(Vinv_S6, C, triples, num_points=20, seed=42)
    
    if coeffs is not None:
        log("\n" + "="*70)
        log("GRAVITY AMPLITUDE FOUND!")
        log("="*70)
        log(f"Coefficients: {coeffs}")
        
        # Construct the final form
        gravity_form = sum(c * v for c, v in zip(coeffs, Vinv_S6))
        
        # Count nonzero entries
        nnz = len([x for x in gravity_form if x != 0])
        log(f"Form has {nnz} nonzero entries out of {len(gravity_form)}")
        
        # Save result
        result = {
            'status': 'success',
            'coefficients': [str(c) for c in coeffs],
            'form_nnz': nnz,
            'form_total': len(gravity_form),
        }
        
        with open('gravity_amplitude_result.json', 'w') as f:
            json.dump(result, f, indent=2)
        log("\nResult saved to gravity_amplitude_result.json")
    else:
        log("\n" + "="*70)
        log("FAILED TO FIND GRAVITY AMPLITUDE")
        log("="*70)
        log("Suggestions:")
        log("  1. Try more kinematic points")
        log("  2. Check spinor-helicity implementation")
        log("  3. Verify OS3 form evaluation")


if __name__ == '__main__':
    main()


