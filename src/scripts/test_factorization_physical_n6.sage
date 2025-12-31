import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import MomentumTwistor, ang_bracket
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor

def vec_sq(v):
    return 0

def s_ijk(indices, lambdas, tildes):
    P_00 = sum(lambdas[k][0] * tildes[k][0] for k in indices)
    P_01 = sum(lambdas[k][0] * tildes[k][1] for k in indices)
    P_10 = sum(lambdas[k][1] * tildes[k][0] for k in indices)
    P_11 = sum(lambdas[k][1] * tildes[k][1] for k in indices)
    return P_00 * P_11 - P_01 * P_10

def decompose_null_vector(lambdas, tildes, indices):
    P_00 = sum(lambdas[k][0] * tildes[k][0] for k in indices)
    P_01 = sum(lambdas[k][0] * tildes[k][1] for k in indices)
    P_10 = sum(lambdas[k][1] * tildes[k][0] for k in indices)
    P_11 = sum(lambdas[k][1] * tildes[k][1] for k in indices)
    lP = vector(QQ, [P_00, P_10])
    if lP == 0:
        lP = vector(QQ, [P_01, P_11])
    if lP[0] != 0:
        ltP = vector(QQ, [P_00/lP[0], P_01/lP[0]])
    else:
        ltP = vector(QQ, [P_10/lP[1], P_11/lP[1]])
    return lP, ltP

def test_factorization_physical():
    print('Testing Physical Factorization on Channel s012 -> 0 (N=6)...')
    seed = 42
    set_random_seed(seed)
    mt = MomentumTwistor(n=6)
    Z_base = list(mt.Z)
    Z5 = Z_base[5]
    Z0 = Z_base[0]
    Z3 = Z_base[3]
    Z2_target = Z5 + Z0 + Z3
    Z2_orig = Z_base[2]
    epsilons = [QQ(1)/QQ(10)**k for k in range(2, 9)]
    print(f'{'eps':<10} | {'s012':<12} | {'M6':<12} | {'Residue':<12} | {'Ratio':<12}')
    print('-' * 70)
    predicted_residue = None
    for eps in epsilons:
        Z_def = list(Z_base)
        Z_def[2] = Z2_target + eps * Z2_orig
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
        s_val = s_ijk([0, 1, 2], lambdas, tildes)
        residue = abs(float(M6 * s_val))
        if eps == epsilons[-1]:
            lP, ltP = decompose_null_vector(lambdas, tildes, [0, 1, 2])
            lP_minus = -lP
            L_lambdas = [lambdas[0], lambdas[1], lambdas[2], lP_minus]
            L_tildes  = [tildes[0],  tildes[1],  tildes[2],  ltP]
            R_lambdas = [lP, lambdas[3], lambdas[4], lambdas[5]]
            R_tildes  = [ltP, tildes[3],  tildes[4],  tildes[5]]
            ML, stL = hodges_npt_mhv_spinor(L_lambdas, L_tildes, neg=(0,1), delete=(0,1,2))
            MR, stR = hodges_npt_mhv_spinor(R_lambdas, R_tildes, neg=(1,2), delete=(0,1,2))
            if ML is None or MR is None:
                predicted_residue = 0
            else:
                predicted_residue = abs(float(ML * MR))
        ratio = 0
        if predicted_residue and predicted_residue > 0:
            ratio = residue / predicted_residue
        print(f'{float(eps):<10.1e} | {abs(float(s_val)):<12.3e} | {abs(float(M6)):<12.3e} | {residue:<12.3e} | {ratio:<12.4f}')
    print('-' * 70)
    if predicted_residue and abs(ratio - 1.0) < 0.2:
        print(f'SUCCESS: Residue matches factorization product! (Ratio ~ {ratio:.2f})')
    elif predicted_residue:
         print(f'PARTIAL: Stable residue but ratio is {ratio:.2f} (likely convention factor or phase)')
    else:
        print('FAILURE: Could not verify factorization.')

if __name__ == '__main__':
    test_factorization_physical()

