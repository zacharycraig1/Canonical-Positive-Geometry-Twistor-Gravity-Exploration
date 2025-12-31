import sys
import os
from sage.all import *

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian, laplacian_from_weights, tree_sum_kirchhoff
from src.chy_oracle.kinematics_samples import sample_twistor, MomentumTwistor

def analyze_poles():
    print("Starting Phase E2: Spinor Valuation Analysis...")
    
    # Analyze poles for <0 1> -> 0
    i, j = 0, 1
    print(f"\nAnalyzing limit <{i} {j}> -> 0")
    
    # 1. Start from generic seed
    seed = 42
    while True:
        twistor = sample_twistor(seed=seed, n=6)
        lambdas_base = [twistor.get_lambda(k) for k in range(6)]
        tildes_base = [twistor.get_tilde_lambda(k) for k in range(6)]
        if all(t is not None for t in tildes_base):
            break
        seed += 1
        
    print(f"Using seed {seed}")
    
    # Reference spinors for weighted Laplacian
    x = vector(QQ, [1, 0])
    y = vector(QQ, [0, 1])
    
    # Perturbation parameter
    # We use numerical fit for slope: log|f| vs log|eps|
    epsilons = [RR(10)**(-k) for k in range(5, 15)]
    
    vals_M = []
    vals_Tree_Plain = []
    vals_Tree_Weighted_3x3 = []
    
    # Vector v not proportional to lambda_i
    # lambda_i is [u, v]. We want <i v_vec> != 0.
    # If lambda_i = [1, 0], use v=[0, 1].
    li = lambdas_base[i]
    v_vec = vector(QQ, [0, 1])
    if li[0]*v_vec[1] - li[1]*v_vec[0] == 0:
        v_vec = vector(QQ, [1, 0])
        
    print(f"  Base lambda_{i}: {li}")
    print(f"  Perturbation direction: {v_vec}")
    
    for eps in epsilons:
        # Deform lambda_j -> lambda_i + eps * v_vec
        # This makes <i j> = eps * <i v> -> 0
        
        lambdas_def = [v for v in lambdas_base] # Copy
        # We want lambda_j to approach lambda_i
        # But wait, if lambda_j -> lambda_i, then <i j> -> 0.
        # Set lambda_j = lambda_i + eps * v_vec
        lambdas_def[j] = lambdas_base[i] + eps * v_vec
        
        # Use FIXED tilde_lambdas (off-shell momentum)
        tildes = tildes_base
        
        # 1. Physical Amplitude
        M, status = hodges_6pt_mhv_spinor(lambdas_def, tildes)
        if status != "ok":
            print(f"  Error at eps={eps}: {status}")
            vals_M.append(0)
        else:
            # Remove helicity factor <0 1>^8 to see A_6 scaling
            h_factor = ang_bracket(lambdas_def[0], lambdas_def[1])**8
            vals_M.append(abs(M / h_factor))
            
        # 2. Plain Tree Sum
        # Recompute weights w_ab = [ab]/<ab>
        weights_plain = {}
        for a in range(6):
            for b in range(a+1, 6):
                ang = ang_bracket(lambdas_def[a], lambdas_def[b])
                sq = tildes[a][0]*tildes[b][1] - tildes[a][1]*tildes[b][0]
                if ang != 0:
                    weights_plain[(a, b)] = sq / ang
                else:
                    weights_plain[(a, b)] = 0 # Should not happen with eps != 0
        
        tree_plain = tree_sum_kirchhoff(weights_plain, 6)
        vals_Tree_Plain.append(abs(tree_plain))
        
        # 3. Weighted Laplacian 3x3 Minor
        L_tilde, C, Phi_tilde = hodges_weighted_laplacian(lambdas_def, tildes, x, y)
        indices_3x3 = [3, 4, 5]
        minor_3x3 = L_tilde.matrix_from_rows_and_columns(indices_3x3, indices_3x3).det()
        vals_Tree_Weighted_3x3.append(abs(minor_3x3))

    # Fit slopes
    import math
    def get_slope(y_vals, x_vals):
        # Linear regression on log-log
        log_x = [math.log(float(x)) for x in x_vals]
        log_y = [math.log(float(y)) for y in y_vals]
        
        # Simple slope between first and last points
        slope = (log_y[-1] - log_y[0]) / (log_x[-1] - log_x[0])
        return slope
    
    slope_M = get_slope(vals_M, epsilons)
    slope_Tree = get_slope(vals_Tree_Plain, epsilons)
    slope_W = get_slope(vals_Tree_Weighted_3x3, epsilons)
    
    print("\nResults:")
    print(f"  Slope M (Amplitude): {slope_M:.4f} (Expected: -2?)")
    print(f"  Slope Tree (Plain):  {slope_Tree:.4f} (Expected: -1)")
    print(f"  Slope Weighted 3x3:  {slope_W:.4f}")
    
    # Check non-adjacent pole <0 2> -> 0
    # ... (similar implementation if needed, but <0 1> is enough to see the gap)
    
    # Non-adjacent check: <0 2> -> 0
    print("\nAnalyzing limit <0 2> -> 0 (Non-adjacent)")
    k = 2
    vals_M_na = []
    vals_Tree_na = []
    
    for eps in epsilons:
        lambdas_def = [v for v in lambdas_base]
        lambdas_def[k] = lambdas_base[i] + eps * v_vec
        
        M, _ = hodges_6pt_mhv_spinor(lambdas_def, tildes_base)
        vals_M_na.append(abs(M))
        
        weights_plain = {}
        for a in range(6):
            for b in range(a+1, 6):
                ang = ang_bracket(lambdas_def[a], lambdas_def[b])
                sq = tildes_base[a][0]*tildes_base[b][1] - tildes_base[a][1]*tildes_base[b][0]
                if ang != 0: weights_plain[(a, b)] = sq / ang
        
        tree_plain = tree_sum_kirchhoff(weights_plain, 6)
        vals_Tree_na.append(abs(tree_plain))

    slope_M_na = get_slope(vals_M_na, epsilons)
    slope_Tree_na = get_slope(vals_Tree_na, epsilons)
    
    # Check weighted minor for non-adjacent
    vals_W_na = []
    for eps in epsilons:
        lambdas_def = [v for v in lambdas_base]
        lambdas_def[k] = lambdas_base[i] + eps * v_vec
        L_tilde, _, _ = hodges_weighted_laplacian(lambdas_def, tildes_base, x, y)
        indices_3x3 = [3, 4, 5]
        vals_W_na.append(abs(L_tilde.matrix_from_rows_and_columns(indices_3x3, indices_3x3).det()))
        
    slope_W_na = get_slope(vals_W_na, epsilons)
    
    print(f"  Slope M (Non-adj):   {slope_M_na:.4f}")
    print(f"  Slope Tree (Non-adj): {slope_Tree_na:.4f}")
    print(f"  Slope Weighted 3x3 (Non-adj): {slope_W_na:.4f}")

    # Check limit <3 4> -> 0 (KEPT indices)
    print("\nAnalyzing limit <3 4> -> 0 (Kept indices)")
    u, v = 3, 4
    vals_M_kept = []
    vals_W_kept = []
    
    for eps in epsilons:
        lambdas_def = [l for l in lambdas_base]
        # Deform lambda_v -> lambda_u
        lu = lambdas_base[u]
        v_vec_k = vector(QQ, [0, 1])
        if lu[0]*v_vec_k[1] - lu[1]*v_vec_k[0] == 0: v_vec_k = vector(QQ, [1, 0])
        
        lambdas_def[v] = lu + eps * v_vec_k
        
        M, _ = hodges_6pt_mhv_spinor(lambdas_def, tildes_base)
        h_f = ang_bracket(lambdas_def[0], lambdas_def[1])**8
        vals_M_kept.append(abs(M/h_f))
        
        L_tilde, _, _ = hodges_weighted_laplacian(lambdas_def, tildes_base, x, y)
        indices_3x3 = [3, 4, 5]
        vals_W_kept.append(abs(L_tilde.matrix_from_rows_and_columns(indices_3x3, indices_3x3).det()))
        
    slope_M_kept = get_slope(vals_M_kept, epsilons)
    slope_W_kept = get_slope(vals_W_kept, epsilons)
    
    print(f"  Slope M (Kept): {slope_M_kept:.4f}")
    print(f"  Slope Weighted 3x3 (Kept): {slope_W_kept:.4f}")

if __name__ == "__main__":
    analyze_poles()

