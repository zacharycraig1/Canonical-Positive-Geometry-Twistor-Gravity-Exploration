
from sage.all import *
from itertools import permutations
import numpy as np

# Load dependencies
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')

# Import diagram logic locally to ensure self-contained script
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

def mandelstam_adapter(twistor, i, j):
    return twistor.get_angle(i, j) * twistor.get_square(i, j)

def s_pole_val(pole, twistor):
    indices = list(pole)
    val = QQ(0)
    for idx1 in range(len(indices)):
        for idx2 in range(idx1+1, len(indices)):
            j = indices[idx1]; k = indices[idx2]
            val += mandelstam_adapter(twistor, j, k)
    return val

def align_basis(seed=42):
    print(f"Attempting to align Bi-adjoint M^-1 with KLT Kernel S (Seed {seed})...")
    
    # 1. Sample Kinematics
    np.random.seed(seed)
    Z_list = []
    for i in range(6):
        z = vector(QQ, [np.random.randint(-50, 51) for _ in range(4)])
        Z_list.append(z)
    twistor = MomentumTwistor(n=6, Z=Z_list, check_domain=False)
    twistor._compute_brackets()
    
    # Check if valid
    for i in range(6):
        if twistor.get_angle(i, (i+1)%6) == 0:
            print("Invalid kinematics (pole)")
            return

    # 2. Compute M (Bi-adjoint)
    permuted_set = [1, 2, 3]
    basis_perms = sorted(list(permutations(permuted_set)))
    
    alpha_diagrams = {}
    beta_diagrams = {}
    for alpha in basis_perms:
        alpha_diagrams[alpha] = get_triangulations([0] + list(alpha) + [4, 5])
    for beta in basis_perms:
        beta_diagrams[beta] = get_triangulations([0] + list(beta) + [5, 4])
        
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
                    s = s_pole_val(pole, twistor)
                    if s == 0: den = 0; break
                    den *= s
                if den != 0: val += 1/den
            row.append(val)
        m_data.append(row)
    
    M = matrix(QQ, m_data)
    if M.is_singular():
        print("M is singular")
        return
    M_inv = M.inverse()
    
    # 3. Compute S (KLT)
    s_data = []
    for alpha in basis_perms:
        row = []
        for beta in basis_perms:
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), twistor, mandelstam_adapter)
            row.append(val)
        s_data.append(row)
    S = matrix(QQ, s_data)
    
    # 4. Compare M_inv and S
    print("\nComparing S and M_inv:")
    print(f"S[0,0]     = {float(S[0,0]):.4e}")
    print(f"M_inv[0,0] = {float(M_inv[0,0]):.4e}")
    
    ratio = S[0,0] / M_inv[0,0]
    print(f"Ratio (0,0) = {float(ratio):.4e}")
    
    # Normalize
    S_norm = S / S[0,0]
    Minv_norm = M_inv / M_inv[0,0]
    
    diff = S_norm - Minv_norm
    max_diff = max([abs(x) for row in diff for x in row])
    print(f"Max Diff (Normalized) = {float(max_diff):.4e}")
    
    if max_diff < 1e-10:
        print("PERFECT ALIGNMENT FOUND! (Up to scalar)")
        return
        
    # 5. Brute Force Permutation
    # Check if rows of S_norm are permutations of rows of Minv_norm
    print("\nChecking row/col permutations...")
    
    # Try to map rows
    mapping = {}
    used_j = set()
    
    for i in range(6):
        row_s = S_norm[i]
        found = False
        for j in range(6):
            if j in used_j: continue
            
            # Check proportionality
            # row_s vs Minv_norm[j]
            # They are normalized to first element, but first element might not be corresponding.
            # Instead check cosine distance or ratio constancy
            
            # Simple check: Sort absolute values and compare?
            # Or just check if row_s = c * row_m?
            
            # We normalized by [0,0].
            # Let's look at raw rows.
            
            # Check correlation
            v1 = vector(RDF, row_s)
            v2 = vector(RDF, Minv_norm[j])
            
            # Check difference
            d = (v1 - v2).norm()
            if d < 1e-5:
                mapping[i] = j
                used_j.add(j)
                found = True
                print(f"S row {i} matches M_inv row {j} (Direct)")
                break
                
            # Check with sign flip
            d = (v1 + v2).norm()
            if d < 1e-5:
                mapping[i] = (j, -1)
                used_j.add(j)
                found = True
                print(f"S row {i} matches M_inv row {j} (flipped)")
                break
                
    if len(mapping) == 6:
        print("Full row mapping found!")
        print(mapping)
    else:
        print(f"Only mapped {len(mapping)} rows.")
        
    # Check if S is just a transpose of M_inv?
    # S_KLT is symmetric usually?
    if not S.is_symmetric():
        print("S is NOT symmetric.")
        
    if not M_inv.is_symmetric():
        print("M_inv is NOT symmetric.") # Should be symmetric as M is sum of 1/products
        
    # Maybe S = M? (Not inverse)
    # Unlikely, dimensions match inverse.
    
if __name__ == "__main__":
    align_basis()





