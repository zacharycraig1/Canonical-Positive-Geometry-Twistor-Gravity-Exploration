
from sage.all import *
from itertools import permutations, combinations
import os

# Load dependencies
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')

class SpinorHelicityAdapter:
    def __init__(self, lambdas, tilde_lambdas):
        self.n = len(lambdas)
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas
    def get_angle(self, i, j):
        return ang_bracket(self.lambdas, i, j)
    def get_square(self, i, j):
        return sq_bracket(self.tilde_lambdas, i, j)

def mandelstam_adapter(twistor, i, j):
    return twistor.get_angle(i, j) * twistor.get_square(i, j)

# Bi-adjoint Scalar Implementation
# ------------------------------

def get_subsets(n):
    """Generate all proper subsets of {0..n-1} representing potential poles."""
    # We only need subsets of size 2 to n-2.
    # Also S and S^c are the same pole. We normalize to contain 0.
    # Actually, simpler to just use subsets and map to Mandelstam.
    # For n=6, sizes 2,3,4.
    pass

def get_planar_poles(order):
    """
    Return list of poles compatible with the ordering (planar).
    A pole is a subset of indices.
    For n=6, poles are contiguous blocks of length 2, 3, 4 (modulo n).
    We normalize poles to be frozensets.
    """
    n = len(order)
    poles = set()
    
    # Extend order for cyclic slicing
    doubled = list(order) + list(order)
    
    for size in range(2, n-1): # 2, 3, ... n-2
        for i in range(n):
            subset = frozenset(doubled[i : i+size])
            # Normalize: if 0 not in subset? No, s_I = s_{I^c}.
            # Let's verify standard identification.
            # Usually we pick the representation with fewer elements, or containing 0?
            # Let's treat them as sets.
            # But s_{0,1} is same as s_{2,3,4,5}.
            # Let's store the canonical representative (e.g. min element is smallest, or size is small)
            
            if len(subset) > n/2:
                 # Take complement
                 full = frozenset(order)
                 subset = full - subset
            elif len(subset) == n/2:
                # If equal size, use lexicographic check to unique-ify
                # e.g. if 0 in subset, keep it?
                if 0 not in subset:
                    full = frozenset(order)
                    subset = full - subset
            
            poles.add(subset)
            
    return list(poles)

def get_triangulations(order):
    """
    Generate all cubic diagrams (maximal sets of compatible poles) for the order.
    For n=6, each diagram has 3 poles.
    We just iterate all combinations of 3 poles and check compatibility.
    Two poles A, B are compatible if A subset B, B subset A, or A disjoint B (subset of complements).
    Since we use canonical representatives (small sets), "disjoint" is the main check.
    Wait, "crossing" is the issue.
    Condition for A, B compatible:
    They are non-crossing in the cyclic order.
    """
    poles = get_planar_poles(order)
    n = len(order)
    num_propagators = n - 3
    
    # Build compatibility graph
    # Two poles are compatible if they don't cross.
    # "Don't cross": In the cyclic order, the vertices of A and B are not interlaced.
    
    def are_compatible(p1, p2):
        # p1, p2 are frozensets of indices
        if p1 == p2: return True
        if p1.isdisjoint(p2): return True # Disjoint
        if p1.issubset(p2) or p2.issubset(p1): return True # Nested
        return False # Intersecting but not nested -> Crossing?
        # Careful: With canonical reps (small sets), "disjoint" handles well.
        # e.g. (0,1) and (1,2) share 1. Are they compatible?
        # (0,1) and (1,2) -> Union (0,1,2). This is a pole?
        # In a cubic graph, poles must be compatible.
        # Actually, simpler: just check if they cross.
        # A set of poles defines a diagram if all pairs are compatible.
        return False 

    # Better compatibility check:
    # Map indices to positions 0..n-1 in 'order'.
    # A pole is a set of positions.
    # It corresponds to a chord.
    # Chords shouldn't cross.
    
    pos_map = {idx: i for i, idx in enumerate(order)}
    
    def to_cyclic_interval(pole):
        # Convert pole (indices) to interval [start, end] in cyclic order
        # Since they are contiguous, this is well-defined.
        # Find start such that next is in pole...
        # Simpler: sort positions. They should be contiguous modulo n.
        positions = sorted([pos_map[x] for x in pole])
        # Check if contiguous
        is_contig = True
        for k in range(len(positions)-1):
            if positions[k+1] != positions[k]+1:
                is_contig = False
                break
        
        if is_contig:
            return (positions[0], positions[-1])
        
        # If not contiguous, it must wrap around.
        # e.g. [0, 1, 5] for n=6.
        # Then the complement is contiguous in the middle.
        # But we generated poles as contiguous blocks.
        # So canonical representatives might have broken contiguity if we took complement?
        # We only took complement if size > n/2.
        # So size <= n/2.
        # A contiguous block of size <= n/2 can always be represented as contiguous without wrap if we rotate?
        # But positions are fixed 0..n-1.
        # e.g. {0, 5} is contiguous mod 6.
        # positions: [0, 5].
        # It represents chord (5,0).
        return (positions, True) # Just return positions
        
    compatible_pairs = set()
    for i in range(len(poles)):
        for j in range(i+1, len(poles)):
            # Check crossing
            # Two sets A, B.
            # Compatible iff A subset B or B subset A or A disjoint B.
            # This works for any sets?
            # Yes, for laminations on a disk, compatible = nested or disjoint.
            # Intersecting non-nested = crossing.
            p1 = poles[i]
            p2 = poles[j]
            if p1.isdisjoint(p2) or p1.issubset(p2) or p2.issubset(p1):
                compatible_pairs.add((i, j))
                
    # Find cliques of size num_propagators (3)
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(poles)))
    G.add_edges_from(list(compatible_pairs))
    
    cliques = [c for c in nx.find_cliques(G) if len(c) == num_propagators]
    
    # Return list of diagrams (each is tuple of poles)
    diagrams = []
    for c in cliques:
        diagrams.append([poles[x] for x in c])
        
    return diagrams

def compute_s_pole(pole, adapter, mandelstam_func):
    # s_I = (sum p_i)^2
    # If size=2: s_ij
    # If size=3: s_ijk
    pole_list = list(pole)
    val = 0
    # Sum s_ab for distinct a,b in pole
    for i in range(len(pole_list)):
        for j in range(i+1, len(pole_list)):
            val += mandelstam_func(adapter, pole_list[i], pole_list[j])
    return val

def check_biadjoint_relation(num_samples=10):
    print(f"Checking Bi-adjoint Relation over {num_samples} samples...")
    
    # Basis: Permutations of {2,3,4} (indices 1,2,3)
    permuted_set = [1, 2, 3]
    basis_perms = sorted(list(permutations(permuted_set)))
    fixed_legs = [0, 4, 5] # 1, 5, 6
    
    for k in range(num_samples):
        # Sample
        try:
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6)
            adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
        except: continue
        
        # 1. Compute m_inverse (Bi-adjoint matrix)
        # m[alpha, beta] = sum_{diagrams compatible with alpha AND beta} 1 / prod(s_poles)
        
        m_data = []
        
        # Precompute diagrams for each alpha
        # Note: diagrams depend on ordering.
        alpha_diagrams = {}
        for alpha in basis_perms:
            full_order = [0] + list(alpha) + [4, 5]
            # Wait, fixed legs are 0, 4, 5. 
            # In KLT code: A(1, alpha, 5, 6) -> 0, alpha, 4, 5.
            # This is the "Left" ordering.
            alpha_diagrams[alpha] = get_triangulations(full_order)
            
        # For beta (Right ordering):
        # In KLT code: A(1, beta, 6, 5) -> 0, beta, 5, 4.
        # This order matters for planarity!
        beta_diagrams = {}
        for beta in basis_perms:
            full_order = [0] + list(beta) + [5, 4] 
            beta_diagrams[beta] = get_triangulations(full_order)
            
        # Debug: Check diagrams for first alpha
        first_alpha = basis_perms[0]
        # print(f"Alpha {first_alpha} diagrams: {len(alpha_diagrams[first_alpha])}")
        
        # Verify diagram count for first alpha
        n_diags = len(alpha_diagrams[first_alpha])
        if n_diags != 14:
            print(f"WARNING: Expected 14 diagrams, found {n_diags}")
            
        for i, alpha in enumerate(basis_perms):
            row = []
            for j, beta in enumerate(basis_perms):
                # Intersect diagrams
                diags_alpha = set([frozenset(d) for d in alpha_diagrams[alpha]])
                diags_beta = set([frozenset(d) for d in beta_diagrams[beta]])
                
                common = diags_alpha.intersection(diags_beta)
                
                val = QQ(0)
                for diag_poles in common:
                    denom = QQ(1)
                    for pole in diag_poles:
                        s_val = compute_s_pole(pole, adapter, mandelstam_adapter)
                        denom *= s_val
                    if denom == 0:
                        # Should not happen for generic kinematics
                        val = QQ(0)
                        break
                    val += 1/denom
                row.append(val)
            m_data.append(row)
            
        M = matrix(QQ, m_data)
        print(f"M[0,0] = {float(M[0,0])}")
        
        # 2. Compute S_KLT
        s_data = []
        for i, alpha in enumerate(basis_perms):
            row = []
            for j, beta in enumerate(basis_perms):
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_adapter)
                row.append(val)
            s_data.append(row)
        S = matrix(QQ, s_data)
        S_sym = (S + S.transpose()) / 2
        
        # 3. Check M * S_sym
        # Should be Identity * scalar?
        
        Product = M * S_sym
        
        # Check if diagonal/identity
        # Normalize by first element
        ref = Product[0,0]
        if ref == 0:
            print("Zero product element 0,0")
            continue
            
        normalized = Product / ref
        
        # Check if Identity
        is_identity = True
        diff_max = 0
        for r in range(6):
            for c in range(6):
                target = 1 if r == c else 0
                diff = abs(normalized[r,c] - target)
                if diff > 1e-10:
                    is_identity = False
                    diff_max = max(diff_max, diff)
                    
        if is_identity:
            print(f"Sample {k}: Match! M * S = c * I")
        else:
            print(f"Sample {k}: Mismatch. Max diff from I: {float(diff_max)}")
            
        # Check signature of M
        # Normalize M to have order 1 elements to avoid precision issues
        norm_factor = abs(M[0,0]) if M[0,0] != 0 else 1e-20
        M_norm = M / norm_factor
        M_float = M_norm.change_ring(RDF)
        
        try:
            evals = M_float.eigenvalues()
            n_pos = sum(1 for e in evals if e > 1e-5)
            n_neg = sum(1 for e in evals if e < -1e-5)
            n_zero = sum(1 for e in evals if abs(e) <= 1e-5)
            print(f"  Signature of M: ({n_pos}, {n_neg}, {n_zero})")
        except: pass

if __name__ == "__main__":
    check_biadjoint_relation(5)

