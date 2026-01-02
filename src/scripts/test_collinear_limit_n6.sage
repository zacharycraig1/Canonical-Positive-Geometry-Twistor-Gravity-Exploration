import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import MomentumTwistor, ang_bracket
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian

def test_collinear_limit():
    print('Testing Collinear Limit 4 || 5 (N=6)...')
    seed = 999
    set_random_seed(seed)
    mt = MomentumTwistor(n=6)
    Z_base = list(mt.Z)
    Z3 = Z_base[3]
    Z5 = Z_base[5]
    alpha = QQ(1)/QQ(2)
    Z4_target = alpha * Z3 + (1-alpha) * Z5
    Z4_dev = Z_base[4]
    epsilons = [QQ(1)/QQ(10)**k for k in range(2, 9)]
    print(f'{'eps':<10} | {'s45':<12} | {'M6':<12} | {'M6*s45':<12} | {'M6*ang45':<12}')
    print('-' * 70)
    for eps in epsilons:
        Z_def = list(Z_base)
        Z_def[4] = Z4_target + eps * Z4_dev
        mt_eps = MomentumTwistor(n=6, Z=Z_def)
        lambdas = [mt_eps.get_lambda(i) for i in range(6)]
        tildes = []
        possible = True
        for i in range(6):
            lt = mt_eps.get_tilde_lambda(i)
            if lt is None: possible = False; break
            tildes.append(lt)
        if not possible: continue
        x = [1, 0]; y = [0, 1]
        M6, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=[0,1,2])
        if M6 is None: continue
        M6_val = abs(float(M6))
        l4, l5 = lambdas[4], lambdas[5]
        lt4, lt5 = tildes[4], tildes[5]
        ang45 = l4[0]*l5[1] - l4[1]*l5[0]
        sq54 = lt5[0]*lt4[1] - lt5[1]*lt4[0]
        s45 = abs(float(ang45 * sq54))
        prod_s = M6_val * s45
        prod_ang = M6_val * abs(float(ang45))
        print(f'{float(eps):<10.1e} | {s45:<12.3e} | {M6_val:<12.3e} | {prod_s:<12.3e} | {prod_ang:<12.3e}')
    print('-' * 70)
    print('Observation: If M6*s45 is stable, it behaves like 1/s pole (YM-squared-like).')
    print('If M6*ang45 is stable, it behaves like 1/angle (YM-like).')
    print('Gravity MHV collinear limits are often less singular or subtle.')

if __name__ == '__main__':
    test_collinear_limit()

