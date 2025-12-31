
from sage.all import *
from itertools import permutations
import os

# Load dependencies
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')
load('src/hodges.sage') # Needed for MomentumTwistor class

# Bi-adjoint Logic (Copied to ensure self-contained/consistent)
def get_planar_poles(order):
    n = len(order)
    poles = set()
    doubled = list(order) + list(order)
    for size in range(2, n-1):
        for i in range(n):
            subset = frozenset(doubled[i : i+size])
            if len(subset) > n/2:
                 full = frozenset(order)
                 subset = full - subset
            elif len(subset) == n/2:
                if 0 not in subset:
                    full = frozenset(order)
                    subset = full - subset
            poles.add(subset)
    return list(poles)

def get_triangulations(order):
    poles = get_planar_poles(order)
    n = len(order)
    num_propagators = n - 3
    compatible_pairs = set()
    for i in range(len(poles)):
        for j in range(i+1, len(poles)):
            p1 = poles[i]; p2 = poles[j]
            if p1.isdisjoint(p2) or p1.issubset(p2) or p2.issubset(p1):
                compatible_pairs.add((i, j))
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(poles)))
    G.add_edges_from(list(compatible_pairs))
    cliques = [c for c in nx.find_cliques(G) if len(c) == num_propagators]
    return [[poles[x] for x in c] for c in cliques]

def compute_s_pole_explicit(pole, lambdas, tilde_lambdas):
    indices = list(pole)
    val = QQ(0)
    for idx1 in range(len(indices)):
        for idx2 in range(idx1+1, len(indices)):
            j = indices[idx1]
            k = indices[idx2]
            ang = lambdas[j][0]*lambdas[k][1] - lambdas[j][1]*lambdas[k][0]
            sq = tilde_lambdas[j][0]*tilde_lambdas[k][1] - tilde_lambdas[j][1]*tilde_lambdas[k][0]
            val += ang * sq
    return val

def mandelstam_adapter(twistor, i, j):
    # Adapter for KLT function
    return twistor.get_angle(i, j) * twistor.get_square(i, j)

def generate_matrices():
    print("Generating M and S matrices for basis alignment...")
    
    # 1. Find a good point
    # Use explicit spinor sampling for M calculation
    # Use MomentumTwistor for KLT calculation (consistency check needed)
    
    import numpy as np
    # Seed 42 usually works
    try:
        lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6, seed=42)
    except:
        print("Sampling failed")
        return

    # Build Momentum Twistor for KLT
    Z_list = []
    x_current = matrix(QQ, 2, 2, 0)
    for i in range(6):
        lam = lambdas[i]
        til = tilde_lambdas[i]
        mu_0 = x_current[0,0]*lam[0] + x_current[0,1]*lam[1]
        mu_1 = x_current[1,0]*lam[0] + x_current[1,1]*lam[1]
        Z_list.append(vector(QQ, [lam[0], lam[1], mu_0, mu_1]))
        p_matrix = matrix(QQ, 2, 2, [lam[0]*til[0], lam[0]*til[1], lam[1]*til[0], lam[1]*til[1]])
        x_current += p_matrix
        
    twistor = MomentumTwistor(n=6, Z=Z_list, check_domain=False)
    twistor._compute_brackets()
    
    # Basis
    permuted_set = [1, 2, 3] 
    basis_perms = sorted(list(permutations(permuted_set)))
    print(f"Basis size: {len(basis_perms)}")
    for i, p in enumerate(basis_perms):
        print(f"  {i}: {p}")

    # 2. Compute M (Bi-adjoint)
    print("Computing M...")
    alpha_diagrams = {}
    beta_diagrams = {}
    for alpha in basis_perms:
        alpha_diagrams[alpha] = get_triangulations([0] + list(alpha) + [4, 5])
    for beta in basis_perms:
        beta_diagrams[beta] = get_triangulations([0] + list(beta) + [5, 4]) # Note 5,4 for Right
        
    m_data = []
    for alpha in basis_perms:
        row = []
        for beta in basis_perms:
            diags_a = set([frozenset(d) for d in alpha_diagrams[alpha]])
            diags_b = set([frozenset(d) for d in beta_diagrams[beta]])
            common = diags_a.intersection(diags_b)
            val = QQ(0)
            for diag in common:
                den = QQ(1)
                for pole in diag:
                    s = compute_s_pole_explicit(pole, lambdas, tilde_lambdas)
                    if s == 0: den = 0; break
                    den *= s
                if den != 0: val += 1/den
            row.append(val)
        m_data.append(row)
    M = matrix(QQ, m_data)
    
    # 3. Compute S (KLT)
    print("Computing S (KLT)...")
    
    def mandelstam_direct(tw, i, j):
        # Ignore twistor arg, use captured spinors
        ang = lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
        sq = tilde_lambdas[i][0]*tilde_lambdas[j][1] - tilde_lambdas[i][1]*tilde_lambdas[j][0]
        return ang * sq

    s_data = []
    for alpha in basis_perms:
        row = []
        for beta in basis_perms:
            # KLT kernel S[alpha|beta]
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), twistor, mandelstam_direct)
            row.append(val)
        s_data.append(row)
    S = matrix(QQ, s_data)
    
    # Save
    save({'M': M, 'S': S, 'basis': basis_perms}, 'klt_search/matrices_MS.sobj')
    print("Saved to klt_search/matrices_MS.sobj")
    
    # Quick check
    P = M * S
    print("M * S (first row):")
    print(P[0])

if __name__ == "__main__":
    generate_matrices()

