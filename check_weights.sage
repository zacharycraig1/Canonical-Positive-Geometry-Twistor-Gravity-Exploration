
import sys
from sage.all import *
load("correct_klt_proof.sage")

def check_weights():
    print("Checking Homogeneity Weights of Hodges Numerator N(Z)...")
    
    # Denominator function
    def get_denom(twistor):
        # (Prod <i,i+1>)^2
        val = QQ(1)
        for i in range(6):
            val *= twistor.get_angle(i, (i+1)%6)
        return val**2
        
    Z = sample_positive_Z_moment_curve(n=6, seed=42)
    tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    H = hodges_6pt_mhv(tw)[0]
    D = get_denom(tw)
    N = H * D
    print(f"Base N: {float(N):.4e}")
    
    # Check weight for each particle
    weights = []
    for i in range(6):
        # Scale Z[i] -> 2*Z[i]
        Z_scaled = list(Z)
        Z_scaled[i] = [2*x for x in Z[i]]
        
        tw_s = MomentumTwistor(n=6, Z=Z_scaled, check_domain=True)
        H_s = hodges_6pt_mhv(tw_s)[0]
        D_s = get_denom(tw_s)
        N_s = H_s * D_s
        
        ratio = N_s / N
        weight = log(ratio, 2)
        weights.append(weight)
        print(f"Particle {i} Weight: {float(weight):.2f}")
        
    total_weight = sum(weights)
    print(f"Total Weight: {float(total_weight):.2f}")
    
    # Check Overall Scaling (Verification)
    Z_all = [[2*x for x in z] for z in Z]
    tw_all = MomentumTwistor(n=6, Z=Z_all, check_domain=True)
    H_all = hodges_6pt_mhv(tw_all)[0]
    D_all = get_denom(tw_all)
    N_all = H_all * D_all
    
    ratio_all = N_all / N
    weight_all = log(ratio_all, 2)
    print(f"Overall Scaling Weight: {float(weight_all):.2f}")
    
    if abs(weight_all - total_weight) < 0.01:
        print("Consistency Check: PASSED")
    else:
        print("Consistency Check: FAILED")

if __name__ == "__main__":
    check_weights()







