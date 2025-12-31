import sys
import os
from sage.all import *

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor, ang_bracket
from src.chy_oracle.laplacian_bridge import weighted_laplacian
from src.chy_oracle.kinematics_samples import sample_twistor

def analyze_n7_poles():
    print("Starting Phase F4: n=7 Pole Valuation Analysis...")
    
    n = 7
    # Use robust seed
    seed = 100
    while True:
        twistor = sample_twistor(seed=seed, n=n)
        lambdas_base = [twistor.get_lambda(k) for k in range(n)]
        tildes_base = [twistor.get_tilde_lambda(k) for k in range(n)]
        if all(t is not None for t in tildes_base):
            break
        seed += 1
    
    print(f"Using seed {seed}")
    
    # Reference spinors
    x = vector(QQ, [1, 2])
    y = vector(QQ, [3, 1])
    
    # Analyze limit <0 1> -> 0 (Adjacent)
    print("\nAnalyzing limit <0 1> -> 0 (Deleted in minor)")
    # We use delete=(0,1,2) for the reconstruction.
    # So <0 1> is one of the poles in the universal prefactor.
    
    epsilons = [RR(10)**(-k) for k in range(5, 12)]
    
    i, j = 0, 1
    # Perturbation
    li = lambdas_base[i]
    v_vec = vector(QQ, [0, 1])
    if li[0]*v_vec[1] - li[1]*v_vec[0] == 0: v_vec = vector(QQ, [1, 0])
    
    vals_M = []
    vals_Minor = []
    
    for eps in epsilons:
        lambdas_def = [v for v in lambdas_base]
        lambdas_def[j] = lambdas_base[i] + eps * v_vec
        
        # 1. Physical Amplitude
        M, status = hodges_npt_mhv_spinor(lambdas_def, tildes_base, neg=(0,1), delete=(0,1,2))
        
        # Strip helicity factor <0 1>^8 to see bare pole
        h_factor = ang_bracket(lambdas_def[0], lambdas_def[1])**8
        vals_M.append(abs(M / h_factor))
        
        # 2. Weighted Laplacian Minor (roots 0,1,2)
        L_tilde, _, _ = weighted_laplacian(n, lambdas_def, tildes_base, x, y)
        indices_keep = [3, 4, 5, 6] # Keep 3,4,5,6
        minor = L_tilde.matrix_from_rows_and_columns(indices_keep, indices_keep).det()
        vals_Minor.append(abs(minor))
        
    import math
    def get_slope(y_vals):
        log_y = [math.log(float(y)) for y in y_vals]
        log_x = [math.log(float(x)) for x in epsilons]
        return (log_y[-1] - log_y[0]) / (log_x[-1] - log_x[0])
        
    slope_M = get_slope(vals_M)
    slope_Minor = get_slope(vals_Minor)
    
    print(f"  Slope M (stripped): {slope_M:.4f} (Expected -2)")
    print(f"  Slope Minor:        {slope_Minor:.4f} (Expected 0 / Finite)")
    
    # Analyze limit <3 4> -> 0 (Kept in minor)
    print("\nAnalyzing limit <3 4> -> 0 (Kept in minor)")
    u, v = 3, 4
    vals_M_kept = []
    vals_Minor_kept = []
    
    for eps in epsilons:
        lambdas_def = [v for v in lambdas_base]
        lu = lambdas_base[u]
        lambdas_def[v] = lu + eps * v_vec # simplified vector
        
        M, _ = hodges_npt_mhv_spinor(lambdas_def, tildes_base, neg=(0,1), delete=(0,1,2))
        h_f = ang_bracket(lambdas_def[0], lambdas_def[1])**8
        vals_M_kept.append(abs(M / h_f))
        
        L_tilde, _, _ = weighted_laplacian(n, lambdas_def, tildes_base, x, y)
        indices_keep = [3, 4, 5, 6]
        minor = L_tilde.matrix_from_rows_and_columns(indices_keep, indices_keep).det()
        vals_Minor_kept.append(abs(minor))
        
    slope_M_kept = get_slope(vals_M_kept)
    slope_Minor_kept = get_slope(vals_Minor_kept)
    
    print(f"  Slope M (Kept):     {slope_M_kept:.4f} (Expected -1)")
    print(f"  Slope Minor (Kept): {slope_Minor_kept:.4f} (Expected -1)")
    
if __name__ == "__main__":
    analyze_n7_poles()

