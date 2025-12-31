
from sage.all import *
from itertools import permutations
import os

load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/hodges.sage')
load('src/klt.sage')

# Bi-adjoint helper
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

def verify_amplitudes(num_samples=5):
    print("Verifying Amplitudes: KLT vs Intersection vs Hodges")
    
    permuted_set = [1, 2, 3] 
    basis_perms = sorted(list(permutations(permuted_set)))
    
    # Precompute diagrams
    alpha_diagrams = {}
    beta_diagrams = {}
    for alpha in basis_perms:
        alpha_diagrams[alpha] = get_triangulations([0] + list(alpha) + [4, 5])
    for beta in basis_perms:
        beta_diagrams[beta] = get_triangulations([0] + list(beta) + [5, 4])
        
    for k in range(num_samples):
        # 1. Sample
        try:
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6, seed=42+k)
        except: continue
        
        # Build twistor for Hodges/KLT
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
        
        def mandelstam_direct(tw, i, j):
            ang = lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
            sq = tilde_lambdas[i][0]*tilde_lambdas[j][1] - tilde_lambdas[i][1]*tilde_lambdas[j][0]
            return ang * sq
            
        def s_pole_func(pole):
            indices = list(pole)
            val = QQ(0)
            for idx1 in range(len(indices)):
                for idx2 in range(idx1+1, len(indices)):
                    j = indices[idx1]; k = indices[idx2]
                    val += mandelstam_direct(None, j, k)
            return val

        # 2. Compute Matrices
        # S_KLT
        s_data = []
        for alpha in basis_perms:
            row = []
            for beta in basis_perms:
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), None, mandelstam_direct)
                row.append(val)
            s_data.append(row)
        S_KLT = matrix(QQ, s_data)
        
        # M_biadjoint
        m_data = []
        m_bad = False
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
                        s = s_pole_func(pole)
                        if s == 0: den = 0; break
                        den *= s
                    if den != 0: val += 1/den
                    else: m_bad = True
                row.append(val)
            m_data.append(row)
        
        if m_bad: continue
        M = matrix(QQ, m_data)
        if M.is_singular(): continue
        S_INT = M.inverse()
        
        # 3. YM Vectors
        try: num = twistor.get_angle(0, 1)**4
        except: continue
        
        vec_L = []
        for alpha in basis_perms:
            order = [0] + list(alpha) + [4, 5]
            den = QQ(1)
            for i in range(6): den *= twistor.get_angle(order[i], order[(i+1)%6])
            vec_L.append(num/den)
        V_L = vector(QQ, vec_L)
        
        vec_R = []
        for beta in basis_perms:
            order = [0] + list(beta) + [5, 4]
            den = QQ(1)
            for i in range(6): den *= twistor.get_angle(order[i], order[(i+1)%6])
            vec_R.append(num/den)
        V_R = vector(QQ, vec_R)
        
        # 4. Amplitudes
        amp_klt = V_L * S_KLT * V_R
        amp_int = V_L * S_INT * V_R
        amp_hodges, _ = hodges_6pt_mhv(twistor)
        
        if amp_hodges != 0:
            r_klt = amp_klt / amp_hodges
            r_int = amp_int / amp_hodges
            print(f"Sample {k}:")
            print(f"  R_KLT (KLT/Hodges) = {float(r_klt):.4e}")
            print(f"  R_INT (Int/Hodges) = {float(r_int):.4e}")
            print(f"  Rel Diff (KLT vs Int) = {float(abs(amp_klt - amp_int)/abs(amp_klt)):.4e}")

if __name__ == "__main__":
    verify_amplitudes(3)

