#!/usr/bin/env sage
from sage.all import *

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')

def parke_taylor_npt_mhv(twistor, order, neg_helicity=(0, 1)):
    """Parke-Taylor for arbitrary n."""
    n = twistor.n
    if len(order) != n:
        return None
    
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        idx_i = order[i]
        idx_j = order[j]
        bracket = twistor.get_angle(idx_i, idx_j)
        if bracket == 0:
            return None
        denom *= bracket
    
    neg_a, neg_b = neg_helicity
    helicity_factor = twistor.get_angle(neg_a, neg_b)
    if helicity_factor == 0:
        return None
    
    if denom == 0:
        return None
    
    return (helicity_factor ** 4) / denom

# Test n=4
print("Debugging n=4...")

Z = []
for i in range(4):
    t = QQ(i + 1) + QQ(1) / QQ(10)
    z = vector(QQ, [QQ(1), t, t*t, t*t*t])
    Z.append(z)

twistor = MomentumTwistor(n=4, Z=Z, check_domain=True)
print(f"Domain OK: {twistor.domain_ok}")

# Build Phi
n = 4
Phi = matrix(QQ, n, n)
x, y = 0, 3

for i in range(n):
    for j in range(n):
        if i != j:
            ij_ang = twistor.get_angle(i, j)
            ij_sq = twistor.get_square(i, j)
            if ij_sq is not None:
                Phi[i, j] = ij_sq / ij_ang

print(f"Phi[0,1] = {Phi[0,1]}")
print(f"Phi[1,2] = {Phi[1,2]}")

# Diagonal
for i in range(n):
    if i == x or i == y:
        if i == y:
            diag_sum = QQ(0)
            for j in range(n):
                if j != i:
                    diag_sum -= Phi[i, j]
            Phi[i, i] = diag_sum
        else:
            Phi[i, i] = QQ(0)
    else:
        ix_ang = twistor.get_angle(i, x)
        iy_ang = twistor.get_angle(i, y)
        diag_sum = QQ(0)
        for j in range(n):
            if j == i:
                continue
            jx_ang = twistor.get_angle(j, x)
            jy_ang = twistor.get_angle(j, y)
            if jx_ang != 0 and jy_ang != 0:
                contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
                diag_sum -= contrib
        Phi[i, i] = diag_sum

print(f"Phi matrix:")
print(Phi)

# Minor
rows_keep = [1, 2, 3]
cols_keep = [1, 2, 3]
Phi_minor = Phi[rows_keep, cols_keep]
print(f"Phi minor:")
print(Phi_minor)
det_minor = Phi_minor.det()
print(f"Det minor: {det_minor}")

# c factor
a12 = twistor.get_angle(1, 2)
a23 = twistor.get_angle(2, 3)
a31 = twistor.get_angle(3, 1)
c123 = QQ(1) / (a12 * a23 * a31)
print(f"c123 = {c123}")

bar_M4 = (-1) * c123 * det_minor
print(f"bar_M4 = {bar_M4}")

# Test KLT
order_alpha = [3, 0, 1, 2]
A_alpha = parke_taylor_npt_mhv(twistor, order_alpha, neg_helicity=(0, 1))
print(f"A[3,0,1,2] = {A_alpha}")

order_beta = [2, 1, 3, 0]
A_beta = parke_taylor_npt_mhv(twistor, order_beta, neg_helicity=(0, 1))
print(f"A[2,1,3,0] = {A_beta}")

s01 = mandelstam_invariant(twistor, 0, 1)
print(f"s_01 = {s01}")

if A_alpha is not None and A_beta is not None and s01 is not None:
    M4_klt = A_alpha * s01 * A_beta
    print(f"M4_KLT = {M4_klt}")
    if bar_M4 != 0:
        ratio = M4_klt / bar_M4
        print(f"Ratio = {ratio}")
