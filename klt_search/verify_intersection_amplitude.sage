
from sage.all import *
from itertools import permutations
import os

# Load dependencies
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/hodges.sage')

# 1. YM MHV Amplitude Vector
def parke_taylor_den(twistor, order):
    """
    Compute denominator of Parke-Taylor factor for specific ordering.
    PTden(alpha) = <a1 a2><a2 a3>...<an a1>
    """
    denom = QQ(1)
    n = len(order)
    for i in range(n):
        idx_i = order[i]
        idx_j = order[(i+1)%n]
        val = twistor.get_angle(idx_i, idx_j)
        if val == 0: return None
        denom *= val
    return denom

def get_ym_mhv_vector(basis_perms, twistor, neg_legs=(0, 1)):
    """
    Compute vector of YM MHV amplitudes for the given basis permutations.
    A_YM(alpha) = <neg1 neg2>^4 / PTden(alpha)
    
    basis_perms: list of permutations of {2,3,4} (indices 1,2,3)
    Full ordering is 1, alpha, 5, 6 (indices 0, alpha, 4, 5)
    """
    neg1, neg2 = neg_legs
    try:
        num = twistor.get_angle(neg1, neg2)**4
    except:
        return None
        
    vec = []
    # Left basis: 1, alpha, 5, 6 (indices 0, alpha, 4, 5)
    fixed_legs = [0, 4, 5] 
    
    for alpha in basis_perms:
        order = [0] + list(alpha) + [4, 5]
        den = parke_taylor_den(twistor, order)
        if den is None: return None
        vec.append(num / den)
        
    return vector(QQ, vec)

def get_ym_mhv_vector_right(basis_perms, twistor, neg_legs=(0, 1)):
    """
    Compute vector of YM MHV amplitudes for the RIGHT basis.
    A_YM(beta) with order 1, beta, 6, 5 (indices 0, beta, 5, 4)
    """
    neg1, neg2 = neg_legs
    try:
        num = twistor.get_angle(neg1, neg2)**4
    except:
        return None
        
    vec = []
    for beta in basis_perms:
        order = [0] + list(beta) + [5, 4]
        den = parke_taylor_den(twistor, order)
        if den is None: return None
        vec.append(num / den)
        
    return vector(QQ, vec)


# 2. Bi-adjoint Logic (reused)
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
            p1 = poles[i]
            p2 = poles[j]
            if p1.isdisjoint(p2) or p1.issubset(p2) or p2.issubset(p1):
                compatible_pairs.add((i, j))
                
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(poles)))
    G.add_edges_from(list(compatible_pairs))
    cliques = [c for c in nx.find_cliques(G) if len(c) == num_propagators]
    diagrams = []
    for c in cliques:
        diagrams.append([poles[x] for x in c])
    return diagrams

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

def verify_hodges_matching(num_samples=20):
    print(f"Verifying Intersection Amplitude vs Full Hodges (n=6, MHV)...")
    
    permuted_set = [1, 2, 3] 
    basis_perms = sorted(list(permutations(permuted_set)))
    
    # Precompute diagrams
    alpha_diagrams = {} 
    beta_diagrams = {}
    for alpha in basis_perms:
        order = [0] + list(alpha) + [4, 5]
        alpha_diagrams[alpha] = get_triangulations(order)
    for beta in basis_perms:
        order = [0] + list(beta) + [5, 4]
        beta_diagrams[beta] = get_triangulations(order)
    
    ratios = []
    
    for k in range(num_samples):
        # 1. Sample Spinors
        try:
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6)
        except: continue

        # 2. Construct Momentum Twistors Z
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
        
        # 3. Compute M (Bi-adjoint Intersection Matrix)
        m_data = []
        m_singular_diag = False
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
                        s_val = compute_s_pole_explicit(pole, lambdas, tilde_lambdas)
                        if s_val == 0: den = 0; break
                        den *= s_val
                    if den == 0: m_singular_diag = True; break
                    val += 1/den
                if m_singular_diag: break
                row.append(val)
            if m_singular_diag: break
            m_data.append(row)
        if m_singular_diag: continue
            
        M = matrix(QQ, m_data)
        try: 
            S_int = M.inverse()
        except: continue
            
        # 4. Compute YM Vectors
        # Helicity: (0, 1) to match Hodges code if it implicitly assumes something?
        # Hodges code uses <i, i+1> denominators. The numerator isn't explicit?
        # Wait, Hodges amplitude function returns the GRAVITY amplitude.
        # hodges_6pt_mhv returns: det'(Phi) / (prod <i,i+1>)^2 * norm_factor?
        # Let's check src/hodges.sage again.
        # M_6^MHV = det'(Phi) / (prod <i,i+1>)^2
        # det'(Phi) has implicit scaling.
        # Let's just use <0 1>^4 for YM numerator.
        
        V_L = get_ym_mhv_vector(basis_perms, twistor, neg_legs=(0, 1))
        V_R = get_ym_mhv_vector_right(basis_perms, twistor, neg_legs=(0, 1))
        
        if V_L is None or V_R is None: continue
        
        # 5. Compute Intersection Gravity Amplitude
        amp_int = V_L * S_int * V_R
        
        # 6. Compute Full Hodges Amplitude
        amp_hodges, reason = hodges_6pt_mhv(twistor)
        if amp_hodges is None: continue
        
        if amp_hodges != 0:
            ratio = amp_int / amp_hodges
            ratios.append(ratio)
            # print(f"Sample {k}: Ratio = {float(ratio):.4e}")

    # Analysis
    if not ratios:
        print("No valid samples.")
        return

    unique_ratios = sorted(list(set(ratios)))
    print(f"\nResults over {len(ratios)} samples:")
    print(f"Unique ratios found: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1:
        print(f"SUCCESS: Constant Ratio = {unique_ratios[0]}")
    else:
        print("FAILURE: Ratio is not constant.")
        print("Sample ratios:")
        for r in unique_ratios[:5]:
            print(f"  {r} (float: {float(r)})")
            
        # Check scaling dimension again if not constant
        # Ratio1/Ratio2
        r1 = ratios[0]
        r2 = ratios[1]
        rat = float(abs(r1/r2))
        print(f"  Rel Diff (sample): {rat}")

if __name__ == "__main__":
    verify_hodges_matching(10)
