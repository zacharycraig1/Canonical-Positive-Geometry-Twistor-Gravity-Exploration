import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import MomentumTwistor, ang_bracket, sq_bracket
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor

def soft_factor_gravity_0(s_idx, lambdas, tildes, n):
    ls = lambdas[s_idx]
    lts = tildes[s_idx]
    xi = vector(QQ, [1, 1])
    if ang_bracket(ls, xi) == 0:
        xi = vector(QQ, [1, 0])
    S = QQ(0)
    ang_xi_s = ang_bracket(xi, ls)
    for a in range(n):
        if a == s_idx: continue
        la = lambdas[a]
        lta = tildes[a]
        sq_sa = sq_bracket(lts, lta)
        ang_sa = ang_bracket(ls, la)
        ang_xi_a = ang_bracket(xi, la)
        term = (sq_sa * ang_xi_a**2) / (ang_sa * ang_xi_s**2)
        S += term
    return S

def test_soft_limit():
    print('Testing Soft Graviton Limit p5 -> 0 (N=6)...')
    s_idx = 5
    n = 6
    set_random_seed(123)
    mt = MomentumTwistor(n=n)
    Z_base = list(mt.Z)
    epsilons = [QQ(1)/QQ(10)**k for k in range(2, 9)]
    print(f'{'eps':<10} | {'M6':<12} | {'S*M5':<12} | {'Ratio':<12}')
    print('-' * 60)
    Z5_orig = Z_base[5]
    for eps in epsilons:
        Z_def = list(Z_base)
        Z_def[5] = eps * Z5_orig
        mt_eps = MomentumTwistor(n=n, Z=Z_def)
        lambdas = [mt_eps.get_lambda(i) for i in range(n)]
        tildes = []
        possible = True
        for i in range(n):
            lt = mt_eps.get_tilde_lambda(i)
            if lt is None: possible = False; break
            tildes.append(lt)
        if not possible: continue
        x = [1, 0]; y = [0, 1]
        M6, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=[0,1,2])
        if M6 is None: continue
        S = soft_factor_gravity_0(s_idx, lambdas, tildes, n)
        l5 = [lambdas[i] for i in range(5)]
        lt5 = [tildes[i] for i in range(5)]
        M5, st5 = hodges_npt_mhv_spinor(l5, lt5, neg=(0,1), delete=(0,1,2))
        if M5 is None: continue
        prediction = S * M5
        ratio = 0
        if prediction != 0:
            ratio = abs(float(M6)) / abs(float(prediction))
        print(f'{float(eps):<10.1e} | {abs(float(M6)):<12.3e} | {abs(float(prediction)):<12.3e} | {ratio:<12.4f}')
    print('-' * 60)
    if abs(ratio - 1.0) < 0.1:
        print(f'SUCCESS: Soft limit satisfied! Ratio ~ {ratio:.2f}')
    else:
        print('FAILURE: Soft limit check failed.')

if __name__ == '__main__':
    test_soft_limit()

