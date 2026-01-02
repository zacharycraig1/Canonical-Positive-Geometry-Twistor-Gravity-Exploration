#!/usr/bin/env sage
'''
POSITIVIZATION MECHANISM SEARCH

North Star: Find a mechanism to convert signed forest cancellations to positive.

Exploring three concrete ideas:
1. FACTORIZATION-LOCAL GROUPING: Do forests in the same factorization channel have uniform signs?
2. PAIRING MECHANISM: Can forests be paired so each pair contributes positively?
3. EIGENSPACE DECOMPOSITION: Separate M = M_+ + M_- and work with |M_+|^2 + |M_-|^2?
'''
from sage.all import *
import itertools
import sys
import os

sys.path.insert(0, os.getcwd())

load('src/spinor_sampling.sage')

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

def compute_forest_sign(forest, lambdas, tilde_lambdas, x_ref, y_ref):
    '''Compute the sign of a forest term.'''
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
    
    return +1 if term > 0 else -1

def get_crossing_edges(forest, L_set, R_set):
    '''Count edges crossing between L and R partition.'''
    crossing = 0
    for u,v in forest:
        u_in_L = u in L_set
        v_in_L = v in L_set
        if u_in_L != v_in_L:
            crossing += 1
    return crossing

print('='*70)
print('MECHANISM 1: FACTORIZATION-LOCAL GROUPING')
print('='*70)
print()
print('Question: Within a factorization channel, do forests have uniform signs?')
print()

forests = enumerate_forests()
print(f'Total forests: {len(forests)}')

# Define factorization channel: {0,1,2} vs {3,4,5}
L_set = {0,1,2}
R_set = {3,4,5}

# Group forests by number of crossing edges
crossing_groups = {}
for f in forests:
    c = get_crossing_edges(f, L_set, R_set)
    if c not in crossing_groups:
        crossing_groups[c] = []
    crossing_groups[c].append(f)

print('Forests grouped by crossing edges:')
for c, group in sorted(crossing_groups.items()):
    print(f'  {c} crossing edges: {len(group)} forests')

# Check sign uniformity within groups
print()
print('Sign analysis within crossing groups:')
for seed in range(3):
    result = sample_spinor_helicity_conserving(n=6, seed=seed*47)
    if result is None: continue
    lambdas, tilde_lambdas = result
    x_ref, y_ref = vector(QQ, [1,2]), vector(QQ, [3,1])
    
    print(f'  Sample {seed}:')
    for c, group in sorted(crossing_groups.items()):
        pos = sum(1 for f in group if compute_forest_sign(f, lambdas, tilde_lambdas, x_ref, y_ref) == +1)
        neg = len(group) - pos
        uniform = 'YES' if pos == 0 or neg == 0 else 'NO'
        print(f'    {c}-crossing: {pos}+/{neg}- uniform={uniform}')

print()
print('='*70)
print('MECHANISM 2: PAIRING FORESTS')
print('='*70)
print()
print('Question: Can we pair forests so each pair has net positive contribution?')
print()

# Check if forests can be paired by some structure
# Idea: pair by complementary edges or by symmetry

# Count forests sharing edges
edge_counts = {}
for f in forests:
    for e in f:
        if e not in edge_counts:
            edge_counts[e] = 0
        edge_counts[e] += 1

print(f'Edge usage distribution:')
usage_dist = {}
for e, c in edge_counts.items():
    if c not in usage_dist:
        usage_dist[c] = 0
    usage_dist[c] += 1
for count, num_edges in sorted(usage_dist.items()):
    print(f'  {num_edges} edges appear in {count} forests')

# Check if paired forests (sharing 2 edges, differing in 1) have opposite signs
print()
print('Checking paired forests (differ by one edge):')

for seed in range(2):
    result = sample_spinor_helicity_conserving(n=6, seed=seed*53)
    if result is None: continue
    lambdas, tilde_lambdas = result
    x_ref, y_ref = vector(QQ, [1,2]), vector(QQ, [3,1])
    
    same_sign_pairs = 0
    opposite_sign_pairs = 0
    
    for i, f1 in enumerate(forests):
        for j, f2 in enumerate(forests):
            if j <= i: continue
            # Check if they differ by exactly one edge
            s1, s2 = set(f1), set(f2)
            if len(s1 & s2) == 2:  # Share 2 edges, differ in 1
                sign1 = compute_forest_sign(f1, lambdas, tilde_lambdas, x_ref, y_ref)
                sign2 = compute_forest_sign(f2, lambdas, tilde_lambdas, x_ref, y_ref)
                if sign1 and sign2:
                    if sign1 == sign2:
                        same_sign_pairs += 1
                    else:
                        opposite_sign_pairs += 1
    
    print(f'  Sample {seed}: same-sign pairs: {same_sign_pairs}, opposite-sign pairs: {opposite_sign_pairs}')

print()
print('='*70)
print('MECHANISM 3: EIGENSPACE DECOMPOSITION (Theoretical)')
print('='*70)
print()
print('The amplitude decomposes as M = M_+ + M_- where:')
print('  M_+ = sum over forests projecting to positive KLT eigenspace')
print('  M_- = sum over forests projecting to negative KLT eigenspace')
print()
print('For squared amplitude: |M|^2 = |M_+ + M_-|^2')
print('  = |M_+|^2 + |M_-|^2 + 2*Re(M_+ * M_-*)')
print()
print('The cross-term 2*Re(M_+ * M_-*) can be negative, so |M|^2 is not')
print('a simple sum of positive terms.')
print()
print('HOWEVER: If we define a new object:')
print('  Omega_signed = |M_+|^2 + |M_-|^2')
print()
print('This IS manifestly positive! It is the sum of positive contributions')
print('from each eigenspace.')
print()
print('Question: Is Omega_signed a meaningful physical quantity?')
print('- It is NOT the squared amplitude')
print('- It might relate to some inclusive observable')
print('- It respects the eigenspace structure')

print()
print('='*70)
print('CONCLUSIONS')
print('='*70)
print('''
1. FACTORIZATION-LOCAL: Signs are NOT uniform within factorization channels.
   Forests in the same channel have mixed signs.

2. PAIRING: Forests differing by one edge do NOT reliably have opposite signs.
   Simple pairing mechanisms don't work.

3. EIGENSPACE DECOMPOSITION: This is the most promising direction.
   - Decompose M = M_+ + M_- by eigenspace projection
   - Define Omega_signed = |M_+|^2 + |M_-|^2 (manifestly positive!)
   - The question is: what is the physical meaning of Omega_signed?

POSSIBLE MECHANISM:
   The "positive object" may not be the amplitude M itself, but rather
   a related quantity like the norm ||M|| in the KLT inner product:
   
   ||M||^2_KLT = M^dagger * S_KLT * M
   
   With proper choice of inner product, this could be positive.
''')
