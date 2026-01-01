
import sys
from sage.all import *
import itertools
from itertools import combinations, combinations_with_replacement

# Load the Hodges implementation
load("correct_klt_proof.sage")

def get_denom_sq(twistor):
    """Compute PT(Z)^2 = (Prod <i,i+1>)^2"""
    val = QQ(1)
    for i in range(6):
        val *= twistor.get_angle(i, (i+1)%6)
    return val**2

def eval_monomial(monomial, twistor):
    """
    Evaluate a basis monomial at a point.
    monomial: list of factors, e.g. [('2b', (0,1)), ('4b', (0,1,2,3))]
    """
    val = QQ(1)
    for type_tag, indices in monomial:
        if type_tag == '2b':
            val *= twistor.get_angle(*indices)
        elif type_tag == '4b':
            val *= twistor.get_four_bracket(*indices)
    return val

def generate_basis_constraints(target_weights):
    """
    Generate all monomials (multiset of 2-brackets and 4-brackets)
    satisfying the target weights.
    target_weights: list of integers [w0, w1, ..., w5]
    """
    basis = []
    n = len(target_weights)
    max_4b = sum(target_weights) // 4
    
    # 4-brackets indices must be sorted tuples of length 4
    all_4b = list(combinations(range(n), 4))
    
    # Helper to solve 2-bracket graph for remaining weights
    def generate_2b_graphs(current_weights):
        # We need to fill remaining weights with 2-brackets (edges)
        # Sum of remaining must be even
        rem_sum = sum(current_weights)
        if rem_sum % 2 != 0: return []
        if rem_sum == 0: return [[]]
        
        # Determine node with first non-zero weight to enforce ordering/avoid dups?
        # Standard recursive graph generation
        
        # Optimization: Process nodes with small non-zero degree first? 
        # Or just fixed order 0..5
        
        # Find first node that needs edges
        first = -1
        for i in range(n):
            if current_weights[i] > 0:
                first = i
                break
        
        if first == -1: return [[]]
        
        # Must connect 'first' to some 'other' > first
        res = []
        
        # Iterate potential neighbors
        for other in range(first + 1, n):
            if current_weights[other] > 0:
                # Add edge (first, other)
                edge = ('2b', (first, other))
                
                new_weights = list(current_weights)
                new_weights[first] -= 1
                new_weights[other] -= 1
                
                # Recurse
                sub_solutions = generate_2b_graphs(new_weights)
                for sol in sub_solutions:
                    res.append([edge] + sol)
                    
        # To handle multigraphs correctly and avoid duplicates like [(0,1), (0,1)] vs [(0,1), (0,1)] order...
        # The recursion above generates ordered lists of edges.
        # But (0,1) then (0,1) is same as (0,1) then (0,1).
        # We need to handle multiplicity explicitly or use a canonical generation.
        
        # Better approach: Iterate edges with multiplicity.
        # But we don't know how many edges.
        
        # Alternative: Recursive fill by node index.
        # "Satisfy node 0": add edges (0, j) until weight[0]=0.
        return generate_2b_by_node(0, current_weights)

    def generate_2b_by_node(node_idx, w):
        if node_idx == n - 1:
            return [[]] if w[node_idx] == 0 else []
            
        if w[node_idx] == 0:
            return generate_2b_by_node(node_idx + 1, w)
            
        # Need to reduce w[node_idx] to 0 by connecting to j > node_idx
        needed = w[node_idx]
        neighbors = range(node_idx + 1, n)
        
        # Distribute 'needed' connections among neighbors
        # Respecting w[neighbor] limits
        
        res = []
        
        # Generator for partitions of 'needed' into bins with capacities
        def partitions(amount, bins, capacities):
            if not bins:
                if amount == 0: yield []
                return
            
            first_cap = capacities[0]
            # Try taking k items for first bin
            # k <= first_cap AND k <= amount
            for k in range(min(amount, first_cap) + 1):
                for tail in partitions(amount - k, bins[1:], capacities[1:]):
                    yield [k] + tail

        caps = [w[j] for j in neighbors]
        
        for p in partitions(needed, neighbors, caps):
            # p is list of edge counts to neighbors
            edges = []
            next_w = list(w)
            next_w[node_idx] = 0 # Satisfied
            
            for i, count in enumerate(p):
                neigh = neighbors[i]
                next_w[neigh] -= count
                for _ in range(count):
                    edges.append(('2b', (node_idx, neigh)))
            
            # Recurse
            subs = generate_2b_by_node(node_idx + 1, next_w)
            for s in subs:
                res.append(edges + s)
                
        return res

    # Main loop over 4-bracket combinations
    # Since m can be up to ~5, and num 4-brackets is 15.
    # combinations_with_replacement(15, m) grows fast.
    # But constraints on weights are tight.
    # Specifically w[2,3,4,5] are small.
    # Prune 4-brackets early.
    
    # Valid 4-brackets: indices from 0..5
    # Filter 4-brackets that are "heavy" on 3 or 4?
    # w[3]=1 => 3 appears in at most 1 4-bracket.
    # w[4]=1 => 4 appears in at most 1 4-bracket.
    # So we can't have multiple 4-brackets containing 3.
    
    valid_4b_indices = []
    for idx in all_4b:
        # Check if compatible with weights (count <= weight)
        if all(1 <= target_weights[x] for x in idx):
            valid_4b_indices.append(idx)
            
    # Iterate m
    for m in range(max_4b + 1):
        # We need to select m 4-brackets.
        # But respecting max counts.
        # Recursive selection of 4-brackets
        
        def select_4b(start_idx, count, current_w):
            if count == 0:
                # Base case: generate 2-brackets
                g2 = generate_2b_graphs(current_w)
                return g2
            
            res_list = []
            for i in range(start_idx, len(valid_4b_indices)):
                b = valid_4b_indices[i]
                # Can we add b?
                possible = True
                next_w = list(current_w)
                for x in b:
                    if next_w[x] >= 1:
                        next_w[x] -= 1
                    else:
                        possible = False
                        break
                
                if possible:
                    # Recurse (allowing replacement? Yes, pass i)
                    # Wait, if we pass i, we allow replacement.
                    # Does 4-bracket squared make sense? Yes.
                    # But w[3]=1 implies 3 cannot appear twice.
                    # So replacement only allowed if weights allow.
                    
                    sub_res = select_4b(i, count - 1, next_w)
                    term = ('4b', b)
                    for s in sub_res:
                        # Prepend term
                        res_list.append([term] + s)
                        
            return res_list

        basis_m = select_4b(0, m, target_weights)
        basis.extend(basis_m)
        
    return basis

# =============================================================================
# SOLVER
# =============================================================================

def solve_exact_numerator():
    print("Enumerating Basis...")
    target = [8, 8, 2, 1, 1, 2]
    basis = generate_basis_constraints(target)
    print(f"Basis size: {len(basis)}")
    
    # Sort basis by <01> power for readability
    def count_01(monomial):
        c = 0
        for type, idx in monomial:
            if type == '2b' and idx == (0,1):
                c += 1
        return c
        
    basis.sort(key=count_01, reverse=True)
    
    print("Top basis elements (by <01> power):")
    for b in basis[:5]:
        print(f"  {b}")
        
    # Collection
    num_samples = len(basis) + 20
    print(f"Collecting {num_samples} samples...")
    
    X_rows = []
    y_vals = []
    
    count = 0
    seed = 8000
    
    while count < num_samples:
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        seed += 1
        tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not tw.domain_ok: continue
        
        H = hodges_6pt_mhv(tw)[0]
        if H is None: continue
        
        D_sq = get_denom_sq(tw)
        N_val = H * D_sq
        
        row = []
        for mon in basis:
            row.append(eval_monomial(mon, tw))
            
        X_rows.append(row)
        y_vals.append(N_val)
        count += 1
        if count % 20 == 0: print(f"  {count}...")
        
    print("Building Matrix...")
    M = matrix(QQ, X_rows)
    Y = vector(QQ, y_vals)
    
    print("Solving...")
    try:
        sol = M.solve_right(Y)
        print("Solution Found!")
        
        # Analyze solution
        terms = []
        for i, c in enumerate(sol):
            if c != 0:
                # Format
                parts = []
                for type, idx in basis[i]:
                    if type == '2b':
                        parts.append(f"<{idx[0]}{idx[1]}>")
                    else:
                        parts.append(f"<{idx[0]}{idx[1]}{idx[2]}{idx[3]}>")
                terms.append((c, "".join(parts)))
                
        terms.sort(key=lambda x: abs(x[0]), reverse=True)
        print("\nNumerator N(Z) =")
        for c, s in terms:
            print(f"  {c} * {s}")
            
    except Exception as e:
        print(f"Failed: {e}")
        # rank check
        print(f"Rank: {M.rank()} / {len(basis)}")

if __name__ == "__main__":
    solve_exact_numerator()







