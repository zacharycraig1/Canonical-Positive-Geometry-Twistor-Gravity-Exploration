#!/usr/bin/env sage
from sage.all import *
import os
import numpy as np
import time

def solve_system():
    print("Step 3: Solving Linear System")
    
    if not os.path.exists("klt_search/basis_metadata.sobj"):
        print("Error: Metadata not found.")
        return
        
    basis_data = load("klt_search/basis_metadata.sobj")
    perms = basis_data['perms']
    monomials = basis_data['monomials']
    total_unknowns = basis_data['total_unknowns']
    
    print(f"Total Unknowns: {total_unknowns}")
    
    # Gather Data
    # List all batch files
    import glob
    batch_files = glob.glob("klt_search/data_batch_*.sobj")
    print(f"Found {len(batch_files)} data batches.")
    
    if len(batch_files) == 0:
        print("No data found. Run 02_generate_data.sage first.")
        return

    # Determine total points
    all_points = []
    for f in batch_files:
        points = load(f)
        all_points.extend(points)
        
    num_points = len(all_points)
    print(f"Total Data Points: {num_points}")
    
    if num_points < total_unknowns:
        print(f"WARNING: Underdetermined system ({num_points} < {total_unknowns}). Need more data.")
    
    # Build Matrix A and vector b
    # A has shape (num_points, total_unknowns)
    # Rows are flattened: for each point, we have terms for each (alpha, beta, monomial)
    
    print("Building Matrix A...")
    # We use numpy for constructing the matrix to manage memory
    
    # Mapping: Index = (i_alpha * len(perms) + i_beta) * len(monomials) + i_mon
    len_p = len(perms)
    len_m = len(monomials)
    
    A = np.zeros((num_points, total_unknowns), dtype=np.float64) # Use float for feasibility
    b = np.zeros(num_points, dtype=np.float64)
    
    start_time = time.time()
    
    for i, p in enumerate(all_points):
        target = float(p['target'])
        mon_vals = [float(x) for x in p['mon_vals']]
        pt_vals = [float(x) for x in p['pt_vals']]
        
        b[i] = target
        
        # Fill row
        # This double loop is the bottleneck if done in pure python.
        # Vectorize:
        # The row is essentially OuterProduct(PT, PT) flattened (size 36)
        # Kronecker product with Monomials (size 220)
        
        # PT_pair_vals = [ PT_alpha * PT_beta for alpha in perms for beta in perms ]
        # Row = [ PT_pair_vals[k] * mon_vals[m] for k in 0..35 for m in 0..219 ]
        
        # Let's do this efficiently
        pt_vec = np.array(pt_vals)
        pt_cross = np.outer(pt_vec, pt_vec).flatten() # Size 36
        mon_vec = np.array(mon_vals) # Size 220
        
        # The full row is Kronecker product: pt_cross (x) mon_vec
        # But verify ordering: 
        # Index = k * len_m + m  --> Matches pt_cross[k] * mon_vec[m]
        # This means we group by PT pair first, then iterate monomials.
        
        full_row = np.kron(pt_cross, mon_vec)
        A[i, :] = full_row
        
        if i % 100 == 0:
            print(f"  Processed {i}/{num_points} rows...", end='\r')
            
    print(f"\nMatrix built in {time.time() - start_time:.2f}s")
    
    # Weighted Least Squares
    print("Applying weights for relative error minimization...")
    weights = np.zeros(num_points)
    for i in range(num_points):
        val = abs(b[i])
        if val > 1e-30:
            weights[i] = 1.0 / val
        else:
            weights[i] = 1.0 # Or 0 to ignore?
            
    # Multiply A and b by weights
    # A_w = W @ A (diagonal W)
    # b_w = W @ b
    
    A_w = A * weights[:, np.newaxis]
    b_w = b * weights
    
    print(f"Matrix Shape: {A.shape}")
    
    # Solve Least Squares
    print("Solving Least Squares (Weighted)...")
    solve_start = time.time()
    
    # Use lstsq on weighted system
    x, residuals, rank, s = np.linalg.lstsq(A_w, b_w, rcond=1e-10)
    
    solve_end = time.time()
    print(f"Solved in {solve_end - solve_start:.2f}s")
    
    # Compute manual residual to be sure (on weighted)
    b_pred_w = A_w @ x
    manual_res = np.sum((b_w - b_pred_w)**2)
    print(f"Weighted Residual Sum (RelErr^2): {manual_res}")
    print(f"Avg Relative Error: {np.sqrt(manual_res / num_points)}")
    print(f"Norm of x: {np.linalg.norm(x)}")
    print(f"Max abs x: {np.max(np.abs(x))}")
    
    print(f"Rank: {rank}/{total_unknowns}")
    
    # Analyze Solution
    # Normalize threshold by max(x)
    max_x = np.max(np.abs(x))
    if max_x == 0:
        print("Solution is all zeros.")
        return

    threshold = max_x * 0.01 # 1% of max
    significant_coeffs = np.where(np.abs(x) > threshold)[0]
    print(f"Significant Coefficients (>1% max): {len(significant_coeffs)}")
    
    # Save solution
    np.save("klt_search/solution_x.npy", x)
    
    # Interpret Top Terms
    print("\nTop 10 Terms:")
    indices = np.argsort(np.abs(x))[::-1][:10]
    for idx in indices:
        val = x[idx]
        
        # Decode index
        # idx = (i_a * len_p + i_b) * len_m + i_mon
        
        i_mon = idx % len_m
        remaining = idx // len_m
        i_b = remaining % len_p
        i_a = remaining // len_p
        
        perm_a = perms[i_a]
        perm_b = perms[i_b]
        mon = monomials[i_mon]
        
        print(f"  {val:.4f} * PT{perm_a} * PT{perm_b} * Monomial{mon}")

if __name__ == "__main__":
    solve_system()

