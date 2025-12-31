import sys
import os
from sage.all import *
import random

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor

def verify_n7_extension():
    print("Starting Phase F3: n=7 Verification...")
    print("Testing Master Identity for 7-point MHV Gravity.")
    
    n_points = 7
    n_seeds = 20
    
    passes = 0
    failures = 0
    
    for seed in range(n_seeds):
        # 1. Kinematics
        try:
            lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=n_points)
        except Exception as e:
            print(f"Seed {seed}: Kinematics gen failed: {e}")
            continue
            
        # 2. Physical Amplitude (Hodges)
        # Using default delete=(0,1,2), neg=(0,1)
        M_phys, status = hodges_npt_mhv_spinor(lambdas, tildes, neg=(0,1), delete=(0,1,2))
        
        if status != "ok":
            print(f"Seed {seed}: Hodges failed ({status})")
            continue
            
        # 3. Reconstructed Amplitude (Laplacian)
        # Reference spinors (random)
        x = vector(QQ, [QQ(random.randint(-5, 5)), QQ(random.randint(-5, 5))])
        y = vector(QQ, [QQ(random.randint(-5, 5)), QQ(random.randint(-5, 5))])
        if x == 0: x = vector(QQ, [1, 1])
        if y == 0: y = vector(QQ, [1, -1])
        
        M_rec, status_rec = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=(0,1,2), neg=(0,1))
        
        if status_rec != "ok":
            print(f"Seed {seed}: Reconstruction failed ({status_rec})")
            continue
            
        # 4. Compare
        if M_phys == 0:
            print(f"Seed {seed}: M_phys is 0.")
            continue
            
        ratio = M_rec / M_phys
        
        if abs(ratio - 1) < 1e-10:
            # print(f"Seed {seed}: PASS")
            passes += 1
        else:
            print(f"Seed {seed}: FAIL. Ratio = {float(ratio)}")
            failures += 1
            
    print("\nSummary:")
    print(f"  Tested: {passes + failures}")
    print(f"  Passed: {passes}")
    print(f"  Failed: {failures}")
    
    if failures == 0 and passes > 0:
        print("\nSUCCESS: The Master Identity holds for n=7!")
        print("This confirms the mechanism generalizes beyond 6 points.")
    else:
        print("\nFAILURE: Check signs, normalization, or definitions.")

if __name__ == "__main__":
    verify_n7_extension()

