import sys
import os
from sage.all import *

# Ensure imports
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import MomentumTwistor, sample_twistor
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_canonical, mandelstam_s

def check_collinear_divergence():
    print("N3: Checking Collinear Divergence (0 || 1)...")
    
    # Base kinematics
    tw_base = sample_twistor(seed=100, n=6)
    Z_base = list(tw_base.Z)
    
    # Deform Z1 -> Z0 + eps * Z_X
    # This makes <0 1> ~ eps
    eps = QQ(1)/1000
    
    Z_X = vector(QQ, [1, 2, 3, 4])
    Z_new = list(Z_base)
    Z_new[1] = Z_base[0] + eps * Z_X
    
    tw_eps = MomentumTwistor(n=6, Z=Z_new)
    
    ls = [tw_eps.get_lambda(i) for i in range(6)]
    lts = []
    for i in range(6):
        lt = tw_eps.get_tilde_lambda(i)
        if lt is None:
            print("Singular tilde")
            return
        lts.append(lt)
    
    # Compute Hodges
    # Use negative helicities 2, 3 to avoid 0,1
    val, status = hodges_npt_mhv_canonical(ls, lts, negative_indices=(2, 3))
    
    if status != "ok":
        print(f"Failed: {status}")
        return
        
    print(f"Hodges (eps={float(eps):.1e}): {val.n():.4e}")
    
    # Compute s_01
    s01 = mandelstam_s(ls, lts, 0, 1)
    print(f"s_01: {s01.n():.4e}")
    
    # Splitting function for Gravity 1/s? 
    # Or 1/<01>^2?
    # YM is 1/sqrt(s)? No, 1/<01><10>?
    
    print(f"Val * s_01: {(val * s01).n():.4e}")
    print(f"Val * s_01^2: {(val * s01**2).n():.4e}")
    print(f"Val * <01>: {(val * tw_eps.get_angle(0,1)).n():.4e}")
    
    # We expect some divergence.
    # If Val * s_01 is constant, it's 1/s pole.
    
    # Let's do a second point to check scaling
    eps2 = eps / 10
    Z_new2 = list(Z_base)
    Z_new2[1] = Z_base[0] + eps2 * Z_X
    tw_eps2 = MomentumTwistor(n=6, Z=Z_new2)
    ls2 = [tw_eps2.get_lambda(i) for i in range(6)]
    lts2 = [tw_eps2.get_tilde_lambda(i) for i in range(6)]
    val2, _ = hodges_npt_mhv_canonical(ls2, lts2, negative_indices=(2, 3))
    s01_2 = mandelstam_s(ls2, lts2, 0, 1)
    
    print(f"Hodges (eps={float(eps2):.1e}): {val2.n():.4e}")
    
    ratio_val = val2 / val
    ratio_eps = eps / eps2 # 10
    
    # If val ~ 1/eps^k
    # ratio_val ~ 10^k
    k = log(ratio_val) / log(ratio_eps)
    print(f"Scaling exponent k (val ~ 1/eps^k): {k.n():.4f}")
    
    # Gravity collinear singularity:
    # 1/s singularity? s ~ eps^2 (if both angle and square scale)?
    # Here Z1 -> Z0 implies lambda1 -> lambda0 (angle -> 0).
    # But mu1 -> mu0 implies tilde1 -> tilde0?
    # Yes. So <01> ~ eps, [01] ~ eps. s ~ eps^2.
    # If 1/s, then 1/eps^2 -> k=2.
    # Let's see.

if __name__ == "__main__":
    check_collinear_divergence()

