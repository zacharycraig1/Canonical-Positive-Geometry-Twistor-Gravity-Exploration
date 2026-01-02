import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_twistor
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_canonical, mandelstam_s
from src.chy_oracle.klt import gravity_6pt_mhv_klt

def check_single_point(seed):
    n = 6
    twistor = sample_twistor(seed=seed, n=n)
    
    lambdas = [twistor.get_lambda(i) for i in range(n)]
    tildes = []
    for i in range(n):
        lt = twistor.get_tilde_lambda(i)
        if lt is None: return None, "singular"
        tildes.append(lt)
        
    # Hodges
    M_hodges, status = hodges_npt_mhv_canonical(lambdas, tildes, negative_indices=(0, 1))
    if status != "ok": return None, f"hodges_{status}"
    
    # KLT
    def mandelstam_func(tw, i, j):
        return mandelstam_s(lambdas, tildes, i, j)
        
    M_klt, status_klt = gravity_6pt_mhv_klt(twistor, mandelstam_func)
    if status_klt != "ok": return None, f"klt_{status_klt}"
    
    if M_klt == 0: return None, "klt_zero"
    
    return M_hodges / M_klt, "ok"

def test_pushforward_vs_klt():
    print("M4.2: End-to-End Pushforward Test (Hodges vs KLT)...")
    
    # Test 1 sample first
    ratio, status = check_single_point(42)
    if status != "ok":
        print(f"Sample 42 failed: {status}")
    else:
        print(f"Sample 42 Ratio: {ratio}")
        
    # Batch Test
    print("\nRunning Batch Test (20 samples)...")
    results = {}
    for i in range(20):
        r, s = check_single_point(i + 100)
        if s == "ok":
            if r not in results: results[r] = 0
            results[r] += 1
        else:
            print(f"  Sample {i+100}: {s}")
            
    print("\nRatio Distribution:")
    for r, count in results.items():
        print(f"  {r}: {count} samples")
        
    # Check for consistency
    if len(results) == 1:
        val = list(results.keys())[0]
        if abs(abs(val) - 1) < 1e-9:
             print("PASS: Consistent agreement (up to sign).")
        else:
             print("WARN: Consistent ratio, but not 1/-1.")
    else:
        print("FAIL: Inconsistent ratios.")

if __name__ == "__main__":
    test_pushforward_vs_klt()




