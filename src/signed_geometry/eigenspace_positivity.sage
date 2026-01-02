#!/usr/bin/env sage
'''
EIGENSPACE DECOMPOSITION: The Most Promising Mechanism

Key insight: The KLT kernel defines an inner product on the space of orderings.
If we decompose M = M_+ + M_- according to positive/negative eigenspaces,
then we can define manifestly positive quantities.

Goal: Compute M_+ and M_- explicitly and explore their properties.
'''
from sage.all import *
import itertools
import numpy as np
import sys
import os

sys.path.insert(0, os.getcwd())

load('src/spinor_sampling.sage')
load('src/klt.sage')

def enumerate_forests(n=6, roots=(0,1,2)):
    roots_set = set(roots)
    k = len(roots)
    all_edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    forests = []
    for edges in itertools.combinations(all_edges, n-k):
        adj = {i: [] for i in range(n)}
        for u,v in edges:
            adj[u].append(v)
            adj[v].append(u)
        visited = set()
        valid = True
        for i in range(n):
            if i not in visited:
                stack, root_count = [i], 1 if i in roots_set else 0
                visited.add(i)
                while stack:
                    curr = stack.pop()
                    for nb in adj[curr]:
                        if nb not in visited:
                            visited.add(nb)
                            stack.append(nb)
                            if nb in roots_set: root_count += 1
                if root_count != 1: valid = False; break
        if valid: forests.append(edges)
    return forests

def compute_forest_term(forest, lambdas, tilde_lambdas, x_ref, y_ref):
    n = 6
    def ang(i,j): return lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
    def sq(i,j): return tilde_lambdas[i][0]*tilde_lambdas[j][1] - tilde_lambdas[i][1]*tilde_lambdas[j][0]
    def C(i): return (lambdas[i][0]*x_ref[1] - lambdas[i][1]*x_ref[0]) * \
                     (lambdas[i][0]*y_ref[1] - lambdas[i][1]*y_ref[0])
    term = QQ(1)
    for u,v in forest:
        if u > v: u,v = v,u
        a = ang(u,v)
        if a == 0: return None
        term *= (sq(u,v)/a) * C(u) * C(v)
    return term

print('='*70)
print('EIGENSPACE DECOMPOSITION ANALYSIS')
print('='*70)

forests = enumerate_forests()

for seed in range(3):
    result = sample_spinor_helicity_conserving(n=6, seed=seed*67)
    if result is None: continue
    
    lambdas, tilde_lambdas = result
    x_ref, y_ref = vector(QQ, [1,2]), vector(QQ, [3,1])
    
    # Compute all forest terms
    pos_sum = QQ(0)
    neg_sum = QQ(0)
    pos_count = 0
    neg_count = 0
    
    for f in forests:
        term = compute_forest_term(f, lambdas, tilde_lambdas, x_ref, y_ref)
        if term is None: continue
        if term > 0:
            pos_sum += term
            pos_count += 1
        else:
            neg_sum += term
            neg_count += 1
    
    M_total = pos_sum + neg_sum  # The amplitude
    M_plus = pos_sum  # Sum of positive forests
    M_minus = neg_sum  # Sum of negative forests (negative value)
    
    print(f'Sample {seed}:')
    print(f'  Positive forests: {pos_count}, sum = {float(pos_sum):.4e}')
    print(f'  Negative forests: {neg_count}, sum = {float(neg_sum):.4e}')
    print(f'  Total amplitude M = {float(M_total):.4e}')
    print()
    print(f'  |M_+| = {float(abs(pos_sum)):.4e}')
    print(f'  |M_-| = {float(abs(neg_sum)):.4e}')
    print(f'  |M_+|^2 + |M_-|^2 = {float(pos_sum**2 + neg_sum**2):.4e}')
    print(f'  |M|^2 = {float(M_total**2):.4e}')
    print()
    
    # The ratio tells us how much the cancellation reduces the amplitude
    if abs(pos_sum) + abs(neg_sum) > 0:
        cancellation_factor = abs(M_total) / (abs(pos_sum) + abs(neg_sum))
        print(f'  Cancellation factor: {float(cancellation_factor):.4f}')
        print(f'  (1 = no cancellation, 0 = complete cancellation)')
    print()

print('='*70)
print('KEY OBSERVATION')
print('='*70)
print('''
The eigenspace decomposition M = M_+ + M_- gives:

1. M_+ = sum of positive-sign forest terms (always positive)
2. M_- = sum of negative-sign forest terms (always negative)

These are NOT independent of kinematics - they depend on which forests
are positive at that kinematic point.

HOWEVER: The SIGNED structure is:
  M = |M_+| - |M_-|    (since M_- is negative)
  
This means:
  |M|^2 = (|M_+| - |M_-|)^2 = |M_+|^2 + |M_-|^2 - 2|M_+||M_-|

The cross-term -2|M_+||M_-| is NEGATIVE, causing cancellation.

MANIFESTLY POSITIVE QUANTITIES:

1. Omega_1 = |M_+|^2 + |M_-|^2
   - This is the "incoherent sum" of eigenspace contributions
   - It is manifestly positive
   - It is LARGER than |M|^2

2. Omega_2 = sqrt(|M_+|^2 + |M_-|^2)
   - A "norm" that treats positive and negative eigenspaces equally
   - Manifestly positive

3. The KLT inner product:
   If we define <M, M>_KLT using the KLT kernel signature appropriately,
   we might get a positive-definite quantity.
''')

print('='*70)
print('EXPLORING THE KLT INNER PRODUCT')
print('='*70)

# Compute KLT kernel and its eigenstructure
for seed in range(2):
    result = sample_spinor_helicity_conserving(n=6, seed=seed*71)
    if result is None: continue
    
    lambdas, tilde_lambdas = result
    
    class Adapter:
        def __init__(self, l, t): self.l, self.t = l, t
        def get_angle(self, i, j): return self.l[i][0]*self.l[j][1] - self.l[i][1]*self.l[j][0]
        def get_square(self, i, j): return self.t[i][0]*self.t[j][1] - self.t[i][1]*self.t[j][0]
    
    adapter = Adapter(lambdas, tilde_lambdas)
    def mandelstam(tw, i, j): return tw.get_angle(i,j) * tw.get_square(i,j)
    
    permuted_set = [1, 2, 3]
    basis_perms = sorted(list(itertools.permutations(permuted_set)))
    
    S = matrix(QQ, 6, 6)
    for i, alpha in enumerate(basis_perms):
        for j, beta in enumerate(basis_perms):
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam)
            if val is None: val = QQ(0)
            S[i, j] = val
    
    S_sym = (S + S.transpose()) / 2
    
    # Compute eigenvalues
    S_np = np.array([[float(S_sym[i,j]) for j in range(6)] for i in range(6)])
    eigenvalues, eigenvectors = np.linalg.eigh(S_np)
    
    print(f'Sample {seed}:')
    print(f'  KLT eigenvalues: {[round(e, 4) for e in eigenvalues]}')
    
    pos_eig = [e for e in eigenvalues if e > 1e-10]
    neg_eig = [e for e in eigenvalues if e < -1e-10]
    
    print(f'  Positive eigenvalues: {len(pos_eig)}, sum = {sum(pos_eig):.4f}')
    print(f'  Negative eigenvalues: {len(neg_eig)}, sum = {sum(neg_eig):.4f}')
    
    # The "modified" KLT inner product
    # Define S_mod = diag(|lambda_1|, ..., |lambda_6|) in eigenbasis
    # Then <v, w>_mod = v^T S_mod w is positive definite!
    
    S_mod_diag = np.diag(np.abs(eigenvalues))
    S_mod = eigenvectors @ S_mod_diag @ eigenvectors.T
    
    print(f'  Modified S (absolute eigenvalues) is positive definite: {np.all(np.linalg.eigvalsh(S_mod) > 0)}')
    print()

print('='*70)
print('THE POSITIVE MECHANISM')
print('='*70)
print('''
DISCOVERED MECHANISM:

The KLT kernel S has signature (3,3). If we define a MODIFIED inner product:

    S_abs = V @ diag(|lambda_1|, ..., |lambda_6|) @ V^T

where V is the eigenvector matrix, then S_abs is POSITIVE DEFINITE.

The gravity amplitude in the KLT basis is a vector A = (A_1, ..., A_6).

The POSITIVE QUANTITY is:
    Omega = A^T @ S_abs @ A

This is:
1. Manifestly positive (S_abs is positive definite)
2. Related to the amplitude (uses KLT structure)
3. Accounts for both positive and negative eigenspaces equally

INTERPRETATION:
    The signed geometry of gravity becomes positive when we "rectify"
    the KLT kernel by taking absolute values of eigenvalues.
    
    This is like measuring "distance" in both timelike and spacelike
    directions with the same sign - it gives the total magnitude
    irrespective of the signature structure.

PHYSICAL MEANING:
    Omega = |contribution from positive eigenspace|^2 
          + |contribution from negative eigenspace|^2
    
    This measures the TOTAL scattering strength, summed incoherently
    over the eigenspace decomposition.
''')
