import sys
import os
import json
import time
from sage.all import *

sys.path.append(os.path.join(os.getcwd(), 'src'))
from chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor

# Load MomentumTwistor class - still needed for kinematics generation
# We can import it from a common place or just load hodges.sage
root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def ang_bracket(la, lb):
    return la[0]*lb[1] - la[1]*lb[0]

def sq_bracket(lta, ltb):
    return lta[0]*ltb[1] - lta[1]*ltb[0]

def get_s_spinor(lambdas, tilde_lambdas, i, j):
    """Computes s_ij = <ij>[ji]"""
    ang = ang_bracket(lambdas[i], lambdas[j])
    sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
    return ang * sq

def get_clearing_denominator(lambdas):
    """
    Computes D_cyclic = (prod <i,i+1>)^2
    """
    n = len(lambdas)
    prod_ang = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        prod_ang *= ang_bracket(lambdas[i], lambdas[j])
    return prod_ang**2

def run_valuations():
    print("Phase C3 - Valuation Engine (Improved)")
    
    # 1. Define Divisors to check
    divisors = []
    # All pairs
    for i in range(6):
        for j in range(i+1, 6):
            divisors.append((i, j))
            
    print(f"Checking {len(divisors)} divisors...")
    
    results = []
    
    # Parameters
    num_seeds = 3
    # Use QQbar? No, keep it QQ for speed/stability if possible.
    # Epsilons should be clean rationals.
    epsilons = [QQ(1)/10**k for k in range(3, 8)] # 1e-3 .. 1e-7
    
    seed_counter = 7000
    
    for div in divisors:
        i_idx, j_idx = div
        is_cyclic = ((j_idx - i_idx) % 6 == 1) or ((i_idx - j_idx) % 6 == 1)
        type_str = "cyclic" if is_cyclic else "non-cyclic"
        if div == (0, 1): type_str += "_special"
        
        print(f"\nDivisor <{i_idx},{j_idx}> ({type_str})")
        
        div_results = []
        
        # We need a robust perturbation strategy.
        # <ij> -> eps.
        # If we just add eps*V to Z[j], then <ij> ~ eps * <i,V>.
        # BUT if <i,V> is small or 0, we have issues.
        # Also, we need to ensure OTHER brackets don't vanish.
        
        for trial in range(num_seeds):
            seed = seed_counter
            seed_counter += 1
            
            # Start with generic configuration
            tw = MomentumTwistor(n=6, seed=seed)
            
            # 1. Force <i,j> = 0 EXACTLY (singular point)
            # Make spinor j proportional to spinor i
            Z_sing = [vector(QQ, z) for z in tw.Z]
            
            # Z_sing[j] should have same spinor as Z_sing[i]
            # Z = (lam, mu). 
            # We keep mu_j same as generic (to avoid full collinearity in twistor space, only spinor collinearity)
            # Actually if lam_i ~ lam_j, then <ij>~0.
            
            # Set lambda_j = lambda_i
            Z_sing[j][0] = Z_sing[i][0]
            Z_sing[j][1] = Z_sing[i][1]
            
            # 2. Perturb Z[j] in direction V
            import random
            random.seed(seed)
            V = vector(QQ, [random.randint(-10,10) for _ in range(4)])
            
            # Ensure <i, V> is not 0
            ang_iV = Z_sing[i][0]*V[1] - Z_sing[i][1]*V[0]
            if ang_iV == 0:
                V[1] += 1 # shift slightly

            # If i and j are adjacent (cyclic), then tilde_lambdas will likely fail 
            # because standard reconstruction divides by <i,j>.
            # We need a robust tilde_lambda reconstruction that doesn't blow up on the boundary
            # OR we accept that M blows up and we measure HOW FAST.
            
            # If tilde_lambda reconstruction fails, use fallback?
            # Or use a different twistor parameterization that keeps <i, i+1> != 0?
            # Wait, we are testing precisely the limit <i, i+1> -> 0.
            # In this limit, tilde_lambda[i] ~ 1/<i, i+1>.
            # So we expect 'tilde_fail' if we use the division formula.
            # BUT we are at eps > 0. The division is valid, just large.
            # So why did it fail? "tilde_recon_fail" happens if denom == 0.
            # At eps=1e-3, <i,j> is not 0.
            # Maybe adjacent brackets <j, k> became 0 accidentally?
            # With random V, unlikely.
            
            # Check if tilde_lambda reconstruction code has explicit checks for 0.
            # It does: "if denom == 0: return None"
            # It shouldn't be exactly 0 unless eps is tiny and we underflow, or we are unlucky.
            # Let's add a check: if denom is small, we proceed but expect large numbers.
            
            # We will use check_domain=False in MomentumTwistor to skip the "return None" checks if possible?
            # No, get_tilde_lambda returns None if denom==0.
            
            scaling_data = []
            
            for eps in epsilons:
                # Z(eps) = Z_sing + eps*V acting on index j
                Z_eps = [v for v in Z_sing]
                Z_eps[j] = Z_eps[j] + eps * V
                
                # Create TW
                # Note: Z_eps entries might be symbolic if eps is symbolic variable?
                # We use rational eps, so it stays in QQ.
                tw_eps = MomentumTwistor(n=6, Z=Z_eps, check_domain=False)
                
                # Check actual angle
                ang_ij = tw_eps.get_angle(i_idx, j_idx)
                
                # Compute M
                lambdas = [tw_eps.get_lambda(k) for k in range(6)]
                tilde_lambdas = [tw_eps.get_tilde_lambda(k) for k in range(6)]
                
                # If tilde_lambda reconstruction fails, it means we hit a bad chart point
                # e.g. <k-1, k> = 0 for some k.
                # If our perturbation kills neighbor brackets, we fail.
                
                if any(x is None for x in tilde_lambdas):
                     scaling_data.append({'eps': float(eps), 'M': None, 'error': 'tilde_fail'})
                     continue
                     
                M_val, reason = hodges_6pt_mhv_spinor(lambdas, tilde_lambdas)
                
                if M_val is None:
                    scaling_data.append({'eps': float(eps), 'M': None, 'error': reason})
                    continue
                    
                D_val = get_clearing_denominator(lambdas)
                N_val = M_val * D_val
                
                scaling_data.append({
                    'eps': float(eps),
                    'ang_ij': float(ang_ij),
                    'M': float(abs(M_val)),
                    'D': float(abs(D_val)),
                    'N': float(abs(N_val))
                })

            # Estimate
            # We want stable integer slopes.
            est = {}
            valid = [d for d in scaling_data if d.get('M') is not None and d.get('ang_ij') != 0]
            
            if len(valid) >= 3: # Need a few points
                # Fit line log(y) = k * log(ang) + C
                import math
                
                # Just take first and last valid?
                p1, p2 = valid[-1], valid[0] # smallest eps, largest eps
                
                try:
                    lx1 = math.log(abs(p1['ang_ij']))
                    lx2 = math.log(abs(p2['ang_ij']))
                    den = lx1 - lx2
                    
                    if abs(den) > 1e-9:
                        for key in ['M', 'D', 'N']:
                            ly1 = math.log(p1[key])
                            ly2 = math.log(p2[key])
                            slope = (ly1 - ly2) / den
                            est[f'k_{key}'] = slope
                except Exception as e:
                    # print(f"Math error: {e}")
                    pass
            
            # Log
            k_M = est.get('k_M')
            k_N = est.get('k_N')
            k_M_str = f"{k_M:.2f}" if k_M is not None else "N/A"
            k_N_str = f"{k_N:.2f}" if k_N is not None else "N/A"
            
            print(f"  Seed {seed}: k(M)={k_M_str}, k(N)={k_N_str}")
            
            div_results.append({
                'seed': seed,
                'estimates': est
            })

        results.append({
            'divisor': div,
            'trials': div_results
        })

    with open('phaseC3_valuations.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nSaved phaseC3_valuations.json")

if __name__ == "__main__":
    run_valuations()
