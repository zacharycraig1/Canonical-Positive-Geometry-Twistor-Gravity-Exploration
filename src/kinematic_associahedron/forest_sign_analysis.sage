# forest_sign_analysis.sage
"""
PROJECTIVE SIGN CONDITION ANALYSIS

Key insight: We don't need all w_ij > 0.
We need the 108 forest terms to have CONSISTENT SIGN.

Each forest term is:
    (-1)^3 × w_{e1} × w_{e2} × w_{e3} × C_i × C_j × C_k × ... 
    = -(w_{e1} × w_{e2} × w_{e3}) × (positive C factors)

If we assume C factors are positive (gauge choice), then:
    sign(term) = -sign(w_{e1} × w_{e2} × w_{e3})

For all 108 terms to have the same sign, we need:
    sign(w_{e1} × w_{e2} × w_{e3}) = same for all forests

This gives us 108 sign conditions. But many might be redundant.

This script:
1. Enumerates all 108 forests for n=6 with roots {0,1,2}
2. Finds the sign conditions on w_ij products
3. Determines the minimal set of independent sign constraints
4. Tests if random kinematics can satisfy these constraints
"""

from sage.all import *
from itertools import combinations, product
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical

print("="*70)
print("FOREST SIGN CONDITION ANALYSIS")
print("="*70)

def enumerate_forests(n, roots):
    """
    Enumerate all spanning forests of K_n with trees rooted at 'roots'.
    
    For n=6 with roots={0,1,2}, each forest has:
    - 3 trees, each rooted at one element of {0,1,2}
    - Each non-root vertex (3,4,5) connects to exactly one other vertex
    - Total: 3 edges
    
    A forest is specified by where each non-root points to.
    """
    non_roots = [v for v in range(n) if v not in roots]
    forests = []
    
    # Each non-root can point to any other vertex
    # But the result must form 3 trees (no cycles in the undirected sense)
    
    # For n=6 with 3 roots and 3 non-roots:
    # Each non-root i ∈ {3,4,5} chooses a parent p_i ∈ {0,1,2,3,4,5} \ {i}
    # Valid if the resulting graph is a forest with 3 components
    
    all_vertices = list(range(n))
    
    # Generate all parent assignments
    for parents in product(all_vertices, repeat=len(non_roots)):
        # parents[i] is the parent of non_roots[i]
        valid = True
        edges = []
        
        for i, nr in enumerate(non_roots):
            p = parents[i]
            if p == nr:  # Can't be own parent
                valid = False
                break
            edges.append((min(nr, p), max(nr, p)))  # Undirected edge
        
        if not valid:
            continue
        
        # Check: is this a valid forest?
        # Use union-find to check connectivity
        parent_uf = {v: v for v in range(n)}
        
        def find(x):
            if parent_uf[x] != x:
                parent_uf[x] = find(parent_uf[x])
            return parent_uf[x]
        
        def union(x, y):
            px, py = find(x), find(y)
            if px == py:
                return False  # Cycle!
            parent_uf[px] = py
            return True
        
        is_forest = True
        for (a, b) in edges:
            if not union(a, b):
                is_forest = False
                break
        
        if not is_forest:
            continue
        
        # Check that each root is in a different component
        root_components = [find(r) for r in roots]
        if len(set(root_components)) != len(roots):
            continue
        
        forests.append(tuple(sorted(edges)))
    
    # Remove duplicates
    forests = list(set(forests))
    return forests


def analyze_sign_patterns(forests):
    """
    Analyze which sign conditions on w_ij are needed for consistent forest signs.
    
    Each forest F has edges {e1, e2, e3}.
    The sign contribution is sign(w_{e1} × w_{e2} × w_{e3}).
    
    For all forests to have the same sign, we need:
    sign(Π w_e for e in F) = constant for all F
    
    This is equivalent to: for any two forests F1, F2,
    sign(Π w_e for e in F1) = sign(Π w_e for e in F2)
    
    Which means: sign(Π w_e for e in F1 △ F2) = +1
    where △ is symmetric difference.
    """
    
    print(f"\n[1] Analyzing {len(forests)} forests...")
    
    # Collect all edges
    all_edges = set()
    for F in forests:
        for e in F:
            all_edges.add(e)
    all_edges = sorted(all_edges)
    
    print(f"    Total distinct edges: {len(all_edges)}")
    print(f"    Edges: {all_edges}")
    
    # Build the "edge matrix" over F_2 (mod 2 arithmetic)
    # Row i corresponds to forest i
    # Column j corresponds to edge j
    # Entry = 1 if edge j is in forest i, else 0
    
    edge_to_idx = {e: i for i, e in enumerate(all_edges)}
    
    M = matrix(GF(2), len(forests), len(all_edges))
    for i, F in enumerate(forests):
        for e in F:
            M[i, edge_to_idx[e]] = 1
    
    print(f"\n[2] Edge matrix dimensions: {M.dimensions()}")
    print(f"    Rank (mod 2): {M.rank()}")
    
    # The sign constraints are: for any two forests F1, F2,
    # the product of w_e over their symmetric difference = positive
    # In log space: sum of log|w_e| signs over symmetric difference = even
    
    # The parity constraints form a vector space over F_2
    # The dimension of this space tells us how many independent sign constraints we have
    
    # Consider differences: (row i) - (row 0) for i = 1, ..., n-1
    # These are the parity vectors
    
    # Create difference matrix properly
    diff_rows = matrix(GF(2), len(forests) - 1, len(all_edges))
    for i in range(1, len(forests)):
        diff_rows[i-1, :] = M[i, :] + M[0, :]  # In GF(2), + is same as -
    
    rank_diff = diff_rows.rank()
    
    print(f"    Rank of difference space: {rank_diff}")
    print(f"    This means: {rank_diff} independent sign constraints")
    
    # Find a basis for the constraint space
    # Each basis vector tells us which edges must have product sign = +1
    
    print(f"\n[3] Finding constraint basis...")
    
    pivots = diff_rows.pivots()
    print(f"    Pivot columns: {pivots}")
    print(f"    Corresponding edges:")
    for p in pivots:
        print(f"      Edge {all_edges[p]}")
    
    return all_edges, M, diff_rows


def test_sign_satisfaction(n_samples=100):
    """
    Test how often random kinematics satisfy the sign constraints.
    """
    
    print(f"\n[4] Testing {n_samples} random kinematics for sign satisfaction...")
    
    forests = enumerate_forests(6, (0, 1, 2))
    print(f"    Enumerated {len(forests)} forests")
    
    all_edges = set()
    for F in forests:
        for e in F:
            all_edges.add(e)
    all_edges = sorted(all_edges)
    
    consistent_count = 0
    sign_distribution = {'+': 0, '-': 0, 'mixed': 0}
    
    for seed in range(1, n_samples + 1):
        try:
            kin = SpinorKinematics.random_rational(n=6, seed=seed)
            lambdas = kin.lambdas
            lambdas_tilde = kin.tilde_lambdas
            
            # Compute w_ij for each edge
            def angle_bracket(i, j):
                return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]
            
            def square_bracket(i, j):
                return lambdas_tilde[i][0] * lambdas_tilde[j][1] - lambdas_tilde[i][1] * lambdas_tilde[j][0]
            
            w = {}
            for (i, j) in all_edges:
                ab = angle_bracket(i, j)
                sb = square_bracket(i, j)
                if ab != 0:
                    w[(i, j)] = sb / ab
                else:
                    w[(i, j)] = None
            
            # Check if any w is undefined
            if any(v is None for v in w.values()):
                continue
            
            # Compute sign of each forest term
            forest_signs = []
            for F in forests:
                product = 1
                for e in F:
                    product *= w[e]
                sign = '+' if product > 0 else '-'
                forest_signs.append(sign)
            
            # Check consistency
            if all(s == '+' for s in forest_signs):
                consistent_count += 1
                sign_distribution['+'] += 1
            elif all(s == '-' for s in forest_signs):
                consistent_count += 1
                sign_distribution['-'] += 1
            else:
                sign_distribution['mixed'] += 1
                
        except Exception as e:
            continue
    
    print(f"\n    Results:")
    print(f"    All positive: {sign_distribution['+']}")
    print(f"    All negative: {sign_distribution['-']}")
    print(f"    Mixed signs: {sign_distribution['mixed']}")
    print(f"    Consistent (all same sign): {consistent_count}/{n_samples}")
    
    return consistent_count, sign_distribution


# Run analysis
print("\n" + "="*70)
print("STEP 1: Enumerate forests")
print("="*70)

forests = enumerate_forests(6, (0, 1, 2))
print(f"Found {len(forests)} forests")

# Show a few examples
print("\nFirst 10 forests (edge sets):")
for i, F in enumerate(forests[:10]):
    print(f"  Forest {i}: {F}")

print("\n" + "="*70)
print("STEP 2: Analyze sign structure")
print("="*70)

all_edges, M, diff_rows = analyze_sign_patterns(forests)

print("\n" + "="*70)
print("STEP 3: Test random kinematics")
print("="*70)

consistent, dist = test_sign_satisfaction(n_samples=200)

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if consistent > 0:
    pct = 100 * consistent / 200
    print(f"""
✅ SIGNIFICANT FINDING:
{consistent}/200 = {pct:.1f}% of random kinematics have CONSISTENT forest signs!

This means the positive region (where amplitude has definite sign) is NOT measure zero.

Next step: Identify the explicit inequalities that define this region.
""")
else:
    print("""
⚠️ No consistent sign patterns found in 200 samples.
The positive region may require special kinematics or may not exist.
""")

