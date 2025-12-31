#!/usr/bin/env sage
from sage.all import *
load('src/hodges.sage')
load('src/klt.sage')

def verify_klt_equivalence():
    print("Verifying KLT Formula vs Hodges Oracle")
    print("======================================")
    
    twistor = MomentumTwistor(n=6, check_domain=False)
    
    # Check 10 points
    for i in range(10):
        # 1. Generate random kinematics
        twistor = MomentumTwistor(n=6, check_domain=False)
        
        # 2. Compute Hodges (Oracle)
        try:
            amp_hodges, reason_hodges = hodges_6pt_mhv(twistor)
            if amp_hodges is None: 
                print(f"Point {i}: Hodges Failed ({reason_hodges})")
                continue
            
            # 3. Compute KLT
            # We need mandelstam function
            amp_klt, reason = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
            
            if reason != "ok":
                print(f"Point {i}: KLT Failed ({reason})")
                continue
                
            # Normalization check
            ang01 = twistor.get_angle(0, 1)
            norm = -(ang01**8)
            amp_hodges_norm = amp_hodges * norm
                
            # 4. Compare
            diff = amp_hodges_norm - amp_klt
            if abs(diff) < 1e-8: # Loose tolerance for now
                print(f"Point {i}: MATCH (Normalized Hodges={amp_hodges_norm:.4f}, KLT={amp_klt:.4f})")
            else:
                ratio = amp_hodges_norm / amp_klt if amp_klt != 0 else 0
                print(f"Point {i}: MISMATCH")
                print(f"  Hodges (raw): {amp_hodges}")
                print(f"  Norm Factor:  {norm}")
                print(f"  Hodges (nrm): {amp_hodges_norm}")
                print(f"  KLT:          {amp_klt}")
                print(f"  Diff:         {diff}")
                print(f"  Ratio:        {ratio}")
                
        except Exception as e:
            print(f"Point {i}: Error {e}")
            import traceback
            traceback.print_exc()
            continue

if __name__ == "__main__":
    verify_klt_equivalence()

