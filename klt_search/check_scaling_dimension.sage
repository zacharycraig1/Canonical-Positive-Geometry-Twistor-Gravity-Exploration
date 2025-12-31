
from sage.all import *
from itertools import permutations
import os

load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/hodges.sage')

def check_scaling_dimension_precise(num_samples=3):
    print("Checking Scaling Dimension of Ratio (M_int / M_hodges)...")
    
    # Use conserved spinor sampling which we know works
    try:
        lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6, seed=42)
    except:
        print("Sampling failed")
        return
        
    # Build twistor from spinors
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
    
    r_base = compute_ratio(twistor)
    if r_base is None:
        print("Base ratio failed: M singular or amp zero")
        return
        
    print(f"Base Ratio = {float(r_base)}")
    
    # Scale Z
    t_val = QQ(2)
    Z_scaled = [z * t_val for z in Z_list]
    twistor_t = MomentumTwistor(n=6, Z=Z_scaled, check_domain=False)
    twistor_t._compute_brackets()
    
    r_scaled = compute_ratio(twistor_t)
    if r_scaled is None:
        print("Scaled ratio failed")
        return
        
    print(f"Scaled Ratio (t=2) = {float(r_scaled)}")
    
    scaling_factor = r_scaled / r_base
    print(f"Scaling Factor = {float(scaling_factor)}")
    
    import math
    if scaling_factor > 0:
        D = math.log(float(scaling_factor), 2)
        print(f"Effective Dimension D = {D}")

def compute_ratio(twistor):
    # print("DEBUG: Computing ratio...")
    permuted_set = [1, 2, 3] 
    basis_perms = sorted(list(permutations(permuted_set)))
    
    alpha_diagrams = {}
    beta_diagrams = {}
    for alpha in basis_perms:
        alpha_diagrams[alpha] = get_triangulations([0] + list(alpha) + [4, 5])
    for beta in basis_perms:
        beta_diagrams[beta] = get_triangulations([0] + list(beta) + [5, 4])
        
    def s_pole_func(pole):
        indices = list(pole)
        val = QQ(0)
        for idx1 in range(len(indices)):
            for idx2 in range(idx1+1, len(indices)):
                j = indices[idx1]
                k = indices[idx2]
                
                ang = twistor.get_angle(j, k)
                if ang == 0: return 0 
                sq = twistor.get_square(j, k)
                if sq is None: return 0
                s_jk = ang * sq
                val += s_jk
        return val

    # 3. M Matrix
    m_data = []
    for alpha in basis_perms:
        row = []
        for beta in basis_perms:
            diags_a = set([frozenset(d) for d in alpha_diagrams[alpha]])
            diags_b = set([frozenset(d) for d in beta_diagrams[beta]])
            common = diags_a.intersection(diags_b)
            # print(f"DEBUG: Common diagrams {len(common)}")
            val = QQ(0)
            for diag in common:
                den = QQ(1)
                for pole in diag:
                    s = s_pole_func(pole)
                    if s == 0: 
                        # print(f"DEBUG: Zero pole {pole}")
                        den = 0; break
                    den *= s
                if den != 0: val += 1/den
                else: 
                     return None
            row.append(val)
        m_data.append(row)
    
    M = matrix(QQ, m_data)
    if M.is_zero() or M.is_singular(): 
        # print("DEBUG: Matrix singular")
        return None
        
    try: S_int = M.inverse()
    except: return None
    
    try:
        num = twistor.get_angle(0, 1)**4
    except: return None
    
    vec_L = []
    for alpha in basis_perms:
        order = [0] + list(alpha) + [4, 5]
        den = QQ(1)
        for i in range(6): 
            ang = twistor.get_angle(order[i], order[(i+1)%6])
            if ang == 0: return None
            den *= ang
        vec_L.append(num/den)
    V_L = vector(QQ, vec_L)
    
    vec_R = []
    for beta in basis_perms:
        order = [0] + list(beta) + [5, 4]
        den = QQ(1)
        for i in range(6): 
            ang = twistor.get_angle(order[i], order[(i+1)%6])
            if ang == 0: return None
            den *= ang
        vec_R.append(num/den)
    V_R = vector(QQ, vec_R)
    
    amp_int = V_L * S_int * V_R
    amp_hodges, _ = hodges_6pt_mhv(twistor)
    
    if amp_hodges == 0 or amp_hodges is None: return None
    return amp_int / amp_hodges

# Helpers
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

if __name__ == "__main__":
    check_scaling_dimension_precise()
