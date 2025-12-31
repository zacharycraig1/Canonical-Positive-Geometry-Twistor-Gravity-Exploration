import sys
import os
import numpy as np
from sage.all import *

sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import MomentumTwistor
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket

def test_factorization_channel():
    print("Testing Factorization Channel (n=6)...")
    
    n = 6
    
    # We want <5 0 1 2> = 0 (Planar channel s01)
    # Search for a good seed
    found_seed = False
    
    for seed in range(100):
        set_random_seed(seed)
        
        Z_vals = []
        for i in range(n):
            z = vector(QQ, [QQ(randint(-10, 10)) for _ in range(4)])
            while z == 0:
                 z = vector(QQ, [QQ(randint(-10, 10)) for _ in range(4)])
            Z_vals.append(z)
            
        # Check orthogonality for all
        bad_seed = False
        for i in range(n):
            if Z_vals[i][0] == 0 or Z_vals[i][1] == 0:
                bad_seed = True
                break
        if bad_seed:
            continue
            
        # Try a quick check at epsilon 1e-4
        Z5 = Z_vals[5]
        Z0 = Z_vals[0]
        Z1 = Z_vals[1]
        
        Z2_singular = Z5 + Z0 + Z1
        Z2_generic = Z_vals[2]
        
        eps = QQ(1)/QQ(10000)
        Z_eps = list(Z_vals)
        Z_eps[2] = Z2_singular + eps * Z2_generic
        
        mt_eps = MomentumTwistor(n=n, Z=Z_eps)
        lambdas = [mt_eps.get_lambda(i) for i in range(n)]
        tildes = []
        possible = True
        for i in range(n):
            lt = mt_eps.get_tilde_lambda(i)
            if lt is None:
                possible = False; break
            tildes.append(lt)
        if not possible: continue
        
        x = [1, 0]; y = [0, 1]
        try:
            M_amp, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=[0,1,2])
        except:
            M_amp = None
            
        if M_amp is None: continue
        
        # Check value
        # <5 0 1 2>
        mat = matrix([Z_eps[5], Z_eps[0], Z_eps[1], Z_eps[2]])
        val_bracket = mat.det()
        
        prod = abs(float(M_amp)) * abs(float(val_bracket))
        prod2 = abs(float(M_amp)) * (abs(float(val_bracket))**2)
        
        if prod > 1e-4 or prod2 > 1e-4:
            print(f"Found interesting seed {seed}: Prod={prod:.2e}, Prod2={prod2:.2e}")
            found_seed = True
            break
            
    if not found_seed:
        print("Could not find a seed with non-vanishing residue.")
        exit(1)
        
    print(f"Running full sweep with seed {seed}...")
    
    print("Checking pole scaling as epsilon -> 0 (Channel <5 0 1 2> -> 0)...")
    print(f"{'epsilon':<10} | {'M_amp':<12} | {'Pole':<12} | {'s01':<12} | {'M*P':<12} | {'Max(z_ij)':<12}")
    print("-" * 80)
    
    epsilons = [QQ(1)/QQ(10)**k for k in range(1, 8)]
    results = []
    
    Z5 = Z_vals[5]; Z0 = Z_vals[0]; Z1 = Z_vals[1]
    Z2_singular = Z5 + Z0 + Z1
    Z2_generic = Z_vals[2]
    
    for eps in epsilons:
        Z_eps = list(Z_vals)
        Z_eps[2] = Z2_singular + eps * Z2_generic
        
        mat = matrix([Z_eps[5], Z_eps[0], Z_eps[1], Z_eps[2]])
        val_bracket = mat.det()
        
        mt_eps = MomentumTwistor(n=n, Z=Z_eps)
        lambdas = [mt_eps.get_lambda(i) for i in range(n)]
        tildes = []
        possible = True
        for i in range(n):
            lt = mt_eps.get_tilde_lambda(i)
            if lt is None: possible = False; break
            tildes.append(lt)
        if not possible:
             print(f"{float(eps):<10.1e} | SKIPPED (Spinors)"); continue
             
        x = [1, 0]; y = [0, 1]
        M_amp, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=[0,1,2])
        if M_amp is None:
             print(f"{float(eps):<10.1e} | SKIPPED ({status})"); continue
             
        try:
             z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
             max_z = max(abs(float(v)) for k, v in z_map.items() if str(k).startswith('z_'))
        except: max_z = -1
        
        M_val = abs(float(M_amp))
        pole_val = abs(float(val_bracket))
        
        # Check s01 = (p0+p1)^2 = <0 1>[1 0]
        s01 = abs(float(ang_bracket(lambdas[0], lambdas[1]) * sq_bracket(tildes[0], tildes[1])))
        
        product = M_val * pole_val
        product2 = M_val * (pole_val**2)
        
        print(f"{float(eps):<10.1e} | {M_val:<12.3e} | {pole_val:<12.3e} | {s01:<12.3e} | {product:<12.3e} | {max_z:<12.3e}")
        results.append((product, product2))
        
    last_p1 = results[-1][0]
    last_p2 = results[-1][1]
    prev_p2 = results[-2][1]
    
    def check_stable(curr, prev):
        if abs(prev) < 1e-15: return False
        return abs((curr - prev)/prev) < 0.1

    print("-" * 75)
    if check_stable(last_p2, prev_p2) and last_p2 > 1e-9:
        print(f"SUCCESS: Amplitude scales as 1/pole^2 (Double Pole). Residue ~ {last_p2:.6e}")
    elif results[-1][0] > 1e-9 and check_stable(results[-1][0], results[-2][0]):
         print(f"SUCCESS: Amplitude scales as 1/pole. Residue ~ {last_p1:.6e}")
    else:
        print("Analysis: Pole seems to vanish or order is unclear.")

if __name__ == "__main__":
    test_factorization_channel()
