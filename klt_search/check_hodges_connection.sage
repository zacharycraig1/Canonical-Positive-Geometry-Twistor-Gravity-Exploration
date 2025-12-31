
from sage.all import *
import os

# Load dependencies
load('src/hodges_sigma.sage')
load('src/hodges.sage')
load('src/klt.sage')

def check_hodges_connection(num_samples=5):
    print(f"Checking Connection between KLT Amplitude and Hodges Amplitude ({num_samples} samples)...")
    
    for k in range(num_samples):
        # 1. Sample Kinematics
        # Use MomentumTwistor which samples Z automatically
        twistor = MomentumTwistor(n=6, check_domain=True)
        if not twistor.domain_ok:
            continue
            
        # 2. Compute KLT Amplitude
        # gravity_6pt_mhv_klt returns (val, status)
        m_klt, status_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        if status_klt != "ok":
            print(f"Sample {k}: KLT failed: {status_klt}")
            continue
            
        # 3. Compute Hodges Amplitude (Full)
        m_hodges, status_hodges = hodges_6pt_mhv(twistor)
        if status_hodges != "ok":
            print(f"Sample {k}: Hodges failed: {status_hodges}")
            continue
            
        # 4. Compute Hodges Reduced Determinant
        m_red, status_red = hodges_6pt_mhv_reduced(twistor)
        
        # 5. Compare
        print(f"\nSample {k}:")
        print(f"  M_KLT    = {m_klt}")
        print(f"  M_Hodges = {m_hodges}")
        print(f"  M_Red    = {m_red}")
        
        if m_hodges != 0:
            ratio = m_klt / m_hodges
            print(f"  Ratio (KLT/Hodges) = {ratio}")
            print(f"  Ratio (float)      = {float(ratio)}")
            
        # Check if Ratio involves <0 1>^8 ?
        # If Ratio is constant, we are good.
        # If Ratio varies, we check scaling.
        
        # Also check scaling of KLT vs Reduced directly
        if m_red != 0:
            ratio_red = m_klt / m_red
            print(f"  Ratio (KLT/Reduced) = {float(ratio_red)}")

    # Check Scaling
    print("\nChecking Scaling Dimension...")
    t = QQ(2)
    twistor = MomentumTwistor(n=6, check_domain=True)
    
    # Base
    klt_1, _ = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
    hodges_1, _ = hodges_6pt_mhv(twistor)
    
    # Scaled
    Z_scaled = [z * t for z in twistor.Z]
    twistor_t = MomentumTwistor(n=6, Z=Z_scaled, check_domain=True)
    klt_t, _ = gravity_6pt_mhv_klt(twistor_t, mandelstam_invariant)
    hodges_t, _ = hodges_6pt_mhv(twistor_t)
    
    print(f"  KLT(t) / KLT(1)       = {float(klt_t/klt_1)} (Expected t^-2 = 0.25)")
    print(f"  Hodges(t) / Hodges(1) = {float(hodges_t/hodges_1)}")
    
    if hodges_1 != 0 and klt_1 != 0:
        ratio_1 = klt_1 / hodges_1
        ratio_t = klt_t / hodges_t
        scaling = ratio_t / ratio_1
        print(f"  Ratio Scaling (t=2)   = {float(scaling)}")
        import math
        dim = math.log(float(scaling), 2)
        print(f"  Ratio Dimension       = {dim}")

if __name__ == "__main__":
    check_hodges_connection(3)

