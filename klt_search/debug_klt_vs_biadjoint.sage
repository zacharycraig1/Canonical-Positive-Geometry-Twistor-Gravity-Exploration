
from sage.all import *
from itertools import permutations
import os

# Load dependencies
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')

# Import helper functions from compute_biadjoint_scalar (re-defined here for standalone)
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

def compute_s_pole(pole, adapter, mandelstam_func):
    pole_list = list(pole)
    val = 0
    for i in range(len(pole_list)):
        for j in range(i+1, len(pole_list)):
            val += mandelstam_func(adapter, pole_list[i], pole_list[j])
    return val

def mandelstam_adapter(twistor, i, j):
    return twistor.get_angle(i, j) * twistor.get_square(i, j)

class SpinorHelicityAdapter:
    def __init__(self, lambdas, tilde_lambdas):
        self.n = len(lambdas)
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas
    def get_angle(self, i, j):
        return ang_bracket(self.lambdas, i, j)
    def get_square(self, i, j):
        return sq_bracket(self.tilde_lambdas, i, j)

def analyze_klt_vs_biadjoint(num_samples=3):
    print(f"Analyzing KLT vs Bi-adjoint Inverse ({num_samples} samples)...")
    
    permuted_set = [1, 2, 3]
    basis_perms = sorted(list(permutations(permuted_set)))
    
    # Cache diagrams
    alpha_diagrams = {}
    beta_diagrams = {}
    for alpha in basis_perms:
        # Left order: 1, alpha, 5, 6 -> indices 0, alpha, 4, 5
        full_order = [0] + list(alpha) + [4, 5]
        alpha_diagrams[alpha] = get_triangulations(full_order)
    for beta in basis_perms:
        # Right order: 1, beta, 6, 5 -> indices 0, beta, 5, 4
        # Note: Standard KLT uses (1, beta, n, n-1)
        full_order = [0] + list(beta) + [5, 4] 
        beta_diagrams[beta] = get_triangulations(full_order)
        
    for k in range(num_samples):
        try:
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6)
            adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
        except: continue
        
        # 1. Compute M (Bi-adjoint)
        m_data = []
        for i, alpha in enumerate(basis_perms):
            row = []
            for j, beta in enumerate(basis_perms):
                diags_alpha = set([frozenset(d) for d in alpha_diagrams[alpha]])
                diags_beta = set([frozenset(d) for d in beta_diagrams[beta]])
                common = diags_alpha.intersection(diags_beta)
                val = QQ(0)
                for diag_poles in common:
                    denom = QQ(1)
                    for pole in diag_poles:
                        s_val = compute_s_pole(pole, adapter, mandelstam_adapter)
                        denom *= s_val
                    val += 1/denom
                row.append(val)
            m_data.append(row)
        M = matrix(QQ, m_data)
        
        # 2. Compute S_KLT
        s_data = []
        for i, alpha in enumerate(basis_perms):
            row = []
            for j, beta in enumerate(basis_perms):
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_adapter)
                row.append(val)
            s_data.append(row)
        S = matrix(QQ, s_data)
        
        # 3. Compute M_inv
        try:
            M_inv = M.inverse()
        except:
            print(f"Sample {k}: M is singular!")
            continue
            
        # 4. Compare S and M_inv
        # Check ratio element-wise
        print(f"\nSample {k}:")
        
        # Check symmetry of S
        print(f"  S symmetric? {S.is_symmetric()}")
        if not S.is_symmetric():
            S_sym = (S + S.transpose())/2
            print(f"  Using S_sym for comparison.")
            S_to_use = S_sym
        else:
            S_to_use = S
            
        # Compare S_to_use and M_inv
        # Ratio of (0,0) elements
        if M_inv[0,0] != 0:
            ratio = S_to_use[0,0] / M_inv[0,0]
            print(f"  Ratio (S_00 / M_inv_00) = {float(ratio)}")
            
            # Check if S ~ ratio * M_inv
            S_pred = ratio * M_inv
            diff = S_to_use - S_pred
            max_diff = max([abs(x) for row in diff for x in row])
            print(f"  Max Diff (S - ratio*M_inv) = {float(max_diff)}")
            
            if max_diff < 1e-10:
                print("  MATCH FOUND! S is proportional to M_inv.")
            else:
                print("  Mismatch.")
                
                # Check eigenvalues/structure
                print("  Diagonal ratio check:")
                for idx in range(6):
                    if M_inv[idx,idx] != 0:
                        r = S_to_use[idx,idx] / M_inv[idx,idx]
                        print(f"    idx {idx}: {float(r)}")
                        
        else:
            print("  M_inv[0,0] is zero.")

if __name__ == "__main__":
    analyze_klt_vs_biadjoint(3)









