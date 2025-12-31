import sys
import os
import random as rnd
import json
import math
from sage.all import *

sys.path.append(os.path.join(os.getcwd(), 'src'))
from chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor
from chy_oracle.bcfw import bcfw_shift_spinors

root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def generate_positive_twistor(n=6, seed=None):
    """
    Generates a positive twistor configuration.
    Positive Geometry means ordered brackets <i i+1 i+2 i+3> > 0.
    Standard construction: Moment curve Z_i(t_i) with ordered t_1 < t_2 < ... < t_n.
    Z(t) = (1, t, t^2, t^3)
    """
    if seed:
        rnd.seed(seed)
        
    # Pick ordered t
    # t_i = i + random
    # Use python range and python ints to avoid Sage types in random.randint
    t_vals = []
    for i in range(1, n+1):
        i_py = int(i)
        t_vals.append(i_py)
    
    Z_list = []
    for i in range(n):
        # Rational approximation
        # Moment curve usually implies R.
        # Let's use rational to stay exact.
        # random.randint(a, b) requires python ints
        start = int(100 * (i+1))
        end = int(100 * (i+2))
        
        val_int = rnd.randint(start, end)
        t_rat = QQ(val_int) / 100
        
        vec = vector(QQ, [1, t_rat, t_rat**2, t_rat**3])
        Z_list.append(vec)
        
    return MomentumTwistor(n=n, Z=Z_list, check_domain=False)

def ang_bracket(la, lb):
    return la[0]*lb[1] - la[1]*lb[0]

def get_clearing_denominator(lambdas):
    n = len(lambdas)
    prod_ang = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        prod_ang *= ang_bracket(lambdas[i], lambdas[j])
    return prod_ang**2

def run_positivity_probe():
    print("Phase C3 - Positivity Probe")
    
    log_data = []
    num_trials = 100
    seed_counter = 10000
    
    pos_N_count = 0
    neg_N_count = 0
    zero_N_count = 0
    
    signs = []
    
    for _ in range(num_trials):
        seed = seed_counter
        seed_counter += 1
        
        # 1. Generate Positive Kinematics
        tw = generate_positive_twistor(n=6, seed=seed)
        
        # 2. Compute Amplitude
        lambdas = [tw.get_lambda(i) for i in range(6)]
        tilde_lambdas = [tw.get_tilde_lambda(i) for i in range(6)]
        
        if any(x is None for x in tilde_lambdas): continue
        
        M_val, reason = hodges_6pt_mhv_spinor(lambdas, tilde_lambdas)
        if M_val is None: continue
        
        # 3. Compute Hat N
        # N = M * D_cyclic
        # Hat N = N / <01>^8 (helicity factor removal)
        
        D_cyclic = get_clearing_denominator(lambdas)
        N_val = M_val * D_cyclic
        
        h_factor = ang_bracket(lambdas[0], lambdas[1])**8
        
        if h_factor == 0: continue
        
        Hat_N = N_val / h_factor
        
        # Check sign
        sign = 0
        if Hat_N > 0: sign = 1
        elif Hat_N < 0: sign = -1
        
        if sign == 1: pos_N_count += 1
        elif sign == -1: neg_N_count += 1
        else: zero_N_count += 1
        
        signs.append(sign)
        
        log_data.append({
            'seed': seed,
            'Hat_N': float(Hat_N),
            'sign': sign
        })
        
    print(f"Positivity Results ({num_trials} trials):")
    print(f"  Positive: {pos_N_count}")
    print(f"  Negative: {neg_N_count}")
    print(f"  Zero: {zero_N_count}")
    
    # Are signs uniform?
    is_uniform = (pos_N_count == 0) or (neg_N_count == 0)
    print(f"  Uniform Sign? {is_uniform}")

    with open('phaseC3_positivity.json', 'w') as f:
        json.dump(log_data, f, indent=2)

if __name__ == "__main__":
    run_positivity_probe()

