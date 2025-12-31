#!/usr/bin/env sage
from sage.all import *
import itertools
import time

# Load libraries
load('src/hodges.sage')  # For MomentumTwistor and hodges_6pt_mhv
load('src/klt.sage')     # For parke_taylor_6pt_mhv and basics

def get_basis_mandelstams(twistor):
    """
    Returns the 9 independent Mandelstam invariants for 6-point.
    Basis: s12, s23, s34, s45, s56, s61, s123, s234, s345
    """
    # 2-particle invariants
    s2 = []
    pairs = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
    for i, j in pairs:
        val = mandelstam_invariant(twistor, i, j)
        if val is None: return None
        s2.append(val)
        
    # 3-particle invariants (s_ijk = s_ij + s_jk + s_ki)
    s3 = []
    triples = [(0,1,2), (1,2,3), (2,3,4)]
    for i, j, k in triples:
        s_ij = mandelstam_invariant(twistor, i, j)
        s_jk = mandelstam_invariant(twistor, j, k)
        s_ki = mandelstam_invariant(twistor, k, i)
        if s_ij is None or s_jk is None or s_ki is None: return None
        s3.append(s_ij + s_jk + s_ki)
        
    return s2 + s3

def generate_monomials(degree, num_vars):
    """
    Generates all monomial exponents for polynomial of given degree.
    Returns list of tuples (e1, ..., en).
    """
    # Using stars and bars / combinations with replacement
    # Simple recursion or itertools
    from itertools import combinations_with_replacement
    indices = range(num_vars)
    monomials = []
    for d in range(degree + 1):
        for combo in combinations_with_replacement(indices, d):
            exponents = [0] * num_vars
            for idx in combo:
                exponents[idx] += 1
            monomials.append(tuple(exponents))
    return monomials

def evaluate_monomial(exponents, values):
    res = 1
    for e, v in zip(exponents, values):
        if e > 0:
            res *= (v ** e)
    return res

def solve_klt_bilinear():
    print("Search for KLT Bilinear Ansatz Coefficients")
    print("===========================================")
    
    # 1. Setup Basis
    # Fixed legs: 1, 5, 6 (indices 0, 4, 5). Permute 2, 3, 4 (indices 1, 2, 3).
    perm_vars = [1, 2, 3]
    perms = sorted(list(itertools.permutations(perm_vars)))
    print(f"Basis size: {len(perms)} x {len(perms)} = {len(perms)**2} terms")
    
    # 2. Setup Polynomial Weights
    # Degree 3 in 9 variables
    degree = 3
    num_vars = 9
    monomial_exps = generate_monomials(degree, num_vars)
    num_coeffs_per_term = len(monomial_exps)
    total_unknowns = (len(perms) ** 2) * num_coeffs_per_term
    
    print(f"Monomials per weight: {num_coeffs_per_term}")
    print(f"Total unknowns: {total_unknowns}")
    
    # We probably want to limit this. 8000 is okay for numerical solve, 
    # but maybe we can restrict to degree 1 first? 
    # KLT Kernel is degree 3 in s (order s^3 for n=6?)
    # Let's check: n=4 -> s^1. n=5 -> s^2. n=6 -> s^3. Yes.
    # So we need degree 3.
    
    # 3. Data Generation
    # We need at least total_unknowns equations.
    num_points = total_unknowns + 50 # Overconstrained
    
    matrix_rows = []
    rhs_vec = []
    
    print(f"Generating {num_points} kinematic points...")
    
    attempts = 0
    collected = 0
    
    while collected < num_points:
        attempts += 1
        if attempts % 100 == 0:
            print(f"  ... collected {collected} points")
            
        # Generate random kinematics
        twistor = MomentumTwistor(n=6, check_domain=False)
        # Re-check domain validity for our specific needs (no division by zero)
        try:
            # Oracle: Use KLT implementation itself to verify ansatz capability
            # We want to rediscover the KLT kernel.
            grav_amp, reason = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
            if reason != "ok": continue
            
            # Basis invariants
            s_vars = get_basis_mandelstams(twistor)
            if s_vars is None: continue
            
            # Monomial values
            mon_vals = [evaluate_monomial(ex, s_vars) for ex in monomial_exps]
            
            # Parke-Taylor values for each basis element
            pt_vals = []
            valid_pt = True
            for p in perms:
                # Ordering: 1, p..., 5, 6 -> [0] + p + [4, 5]
                order = [0] + list(p) + [4, 5]
                val = parke_taylor_6pt_mhv(twistor, order)
                if val is None: 
                    valid_pt = False
                    break
                pt_vals.append(val)
            
            if not valid_pt: continue
            
            # Build row
            # Ansatz = sum_{a,b} (sum_k c_{abk} M_k(s)) * PT_a * PT_b
            #        = sum_{a,b,k} c_{abk} * (M_k(s) * PT_a * PT_b)
            # The row is the vector of coefficients of c_{abk}.
            
            row = []
            for i_a in range(len(perms)):
                for i_b in range(len(perms)):
                    # Term C_ab
                    prefactor = pt_vals[i_a] * pt_vals[i_b]
                    # Expand with monomials
                    term_row = [prefactor * mv for mv in mon_vals]
                    row.extend(term_row)
            
            matrix_rows.append(row)
            rhs_vec.append(grav_amp)
            collected += 1
            
        except Exception as e:
            # print(e)
            continue

    print("Building matrix...")
    # Use simple numerical linear algebra or finite fields if we want exactness
    # For speed and "existence", floating point might be safer first, 
    # but Sage QQ is exact. 8000x8000 in QQ might be slow.
    # Let's try finite field GF(10007) for quick rank check?
    # Or just try to solve a smaller subset first.
    
    # Actually, 8000x8000 dense matrix in QQ is VERY slow.
    # We should use numerical values (RDF) for feasibility check.
    
    # Let's convert to RDF (Real Double Field)
    import numpy as np
    
    print("Converting to numpy/float for speed...")
    # Create numpy arrays
    A = np.array(matrix_rows, dtype=float)
    b = np.array(rhs_vec, dtype=float)
    
    print(f"Matrix shape: {A.shape}")
    
    # Least squares
    print("Solving least squares...")
    start_time = time.time()
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=1e-10)
    end_time = time.time()
    
    print(f"Solved in {end_time - start_time:.2f} seconds")
    print(f"Residuals: {residuals}")
    print(f"Rank: {rank}")
    
    if len(residuals) > 0 and residuals[0] < 1e-10:
        print("SUCCESS: Solution found within numerical precision!")
        print("There exists a bilinear form matching the amplitude.")
        
        # Analyze the solution 'x' to see sparsity or structure
        # x is a flat vector of length total_unknowns
        # Map back to c_{abk}
        
        # Identify non-zero coefficients (threshold)
        threshold = 1e-6
        non_zeros = np.where(np.abs(x) > threshold)[0]
        print(f"Number of non-zero coefficients: {len(non_zeros)}")
        
        # Try to print the top terms
        print("Top 10 terms:")
        indices = np.argsort(np.abs(x))[::-1][:10]
        for idx in indices:
            # Decode idx -> (a, b, k)
            k = idx % num_coeffs_per_term
            ab = idx // num_coeffs_per_term
            b = ab % len(perms)
            a = ab // len(perms)
            
            val = x[idx]
            mon = monomial_exps[k]
            perm_a = perms[a]
            perm_b = perms[b]
            print(f"  Coeff {val:.4f} * PT({perm_a}) * PT({perm_b}) * Monomial{mon}")
            
    else:
        print("FAILURE: No exact solution found.")
        
    return

if __name__ == "__main__":
    solve_klt_bilinear()

