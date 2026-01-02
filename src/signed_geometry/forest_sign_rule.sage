#!/usr/bin/env sage
"""
Forest Sign Rule Discovery

This module investigates what combinatorial property determines the sign
ε(F) = ±1 of each forest term in the gravity amplitude expansion.

The 108 forest terms split 54/54 between positive and negative.
This is not random - there must be a rule.

Candidate hypotheses:
1. Parity of number of edges in certain positions
2. Sign of some determinant associated with the forest
3. Product of edge signs under some orientation
4. Relationship to KLT kernel structure
5. Euler characteristic or other topological invariant
"""
from sage.all import *
import itertools
import sys
import os

sys.path.insert(0, os.getcwd())

load('src/spinor_sampling.sage')


def ang_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def sq_bracket(tilde_lambdas, i, j):
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def enumerate_rooted_forests(n, roots):
    """
    Enumerate all spanning forests of K_n with k trees,
    where each tree contains exactly one root from 'roots'.
    """
    num_roots = len(roots)
    num_edges = n - num_roots
    roots_set = set(roots)
    
    all_edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    
    for edges in itertools.combinations(all_edges, num_edges):
        adj = {i: [] for i in range(n)}
        for u, v in edges:
            adj[u].append(v)
            adj[v].append(u)
        
        visited = set()
        components = []
        valid = True
        
        for i in range(n):
            if i not in visited:
                stack = [i]
                visited.add(i)
                comp = [i]
                root_count = 1 if i in roots_set else 0
                
                while stack:
                    curr = stack.pop()
                    for neighbor in adj[curr]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)
                            comp.append(neighbor)
                            if neighbor in roots_set:
                                root_count += 1
                
                components.append(comp)
                
                if root_count != 1:
                    valid = False
                    break
        
        if valid and len(components) == num_roots:
            yield edges


def compute_forest_sign(forest, lambdas, tilde_lambdas, x_spinor, y_spinor):
    """
    Compute the actual sign of a forest's contribution.
    
    Returns:
        (sign, weight) where sign is +1 or -1
    """
    n = len(lambdas)
    
    def ang_with_ref(lam, ref):
        return lam[0] * ref[1] - lam[1] * ref[0]
    
    # Compute C_i = <i,x><i,y>
    C = [ang_with_ref(lambdas[i], x_spinor) * ang_with_ref(lambdas[i], y_spinor) 
         for i in range(n)]
    
    # Compute weight = prod_{(i,j) in forest} (-w_ij * C_i * C_j)
    weight = QQ(1)
    for u, v in forest:
        ang = ang_bracket(lambdas, u, v)
        sq = sq_bracket(tilde_lambdas, u, v)
        
        if ang == 0:
            return None, None
        
        w_uv = sq / ang
        weight *= (-w_uv * C[u] * C[v])
    
    sign = 1 if weight > 0 else -1 if weight < 0 else 0
    return sign, weight


# ============================================================
# HYPOTHESIS 1: Parity of edges crossing a partition
# ============================================================

def count_crossing_edges(forest, partition):
    """
    Count edges that cross between the two sets of the partition.
    """
    left, right = partition
    left_set = set(left)
    right_set = set(right)
    
    count = 0
    for u, v in forest:
        if (u in left_set and v in right_set) or (u in right_set and v in left_set):
            count += 1
    return count


def hypothesis_crossing_parity(forest, partition):
    """
    Hypothesis: sign = (-1)^{number of crossing edges}
    """
    n_cross = count_crossing_edges(forest, partition)
    return (-1)**n_cross


# ============================================================
# HYPOTHESIS 2: Edge orientation sign
# ============================================================

def hypothesis_edge_orientation(forest):
    """
    Hypothesis: sign = product of (-1)^{i>j} for each edge (i,j)
    
    This counts the number of edges oriented "backwards" in the
    canonical ordering.
    """
    sign = 1
    for u, v in forest:
        if u > v:
            sign *= -1
    return sign


# ============================================================
# HYPOTHESIS 3: Forest structure - inversions
# ============================================================

def compute_forest_inversions(forest, n, roots):
    """
    Compute the number of inversions in the forest structure.
    
    An inversion occurs when a smaller-indexed vertex is a descendant
    of a larger-indexed vertex (relative to the roots).
    """
    roots_set = set(roots)
    
    # Build adjacency
    adj = {i: [] for i in range(n)}
    for u, v in forest:
        adj[u].append(v)
        adj[v].append(u)
    
    inversions = 0
    
    # For each tree, root it and count inversions
    for root in roots:
        # BFS to find parent-child relationships
        parent = {root: None}
        queue = [root]
        
        while queue:
            curr = queue.pop(0)
            for neighbor in adj[curr]:
                if neighbor not in parent:
                    parent[neighbor] = curr
                    queue.append(neighbor)
                    # Count inversion if child < parent
                    if neighbor < curr:
                        inversions += 1
    
    return inversions


def hypothesis_inversion_parity(forest, n, roots):
    """
    Hypothesis: sign = (-1)^{number of inversions}
    """
    inv = compute_forest_inversions(forest, n, roots)
    return (-1)**inv


# ============================================================
# HYPOTHESIS 4: Determinant sign of incidence matrix
# ============================================================

def compute_incidence_matrix(forest, n, roots):
    """
    Build the incidence matrix of the forest.
    
    Rows = non-root vertices
    Cols = edges (oriented from lower to higher index)
    Entry = +1 if edge enters vertex, -1 if leaves, 0 otherwise
    """
    non_roots = [i for i in range(n) if i not in roots]
    edges = list(forest)
    
    m = len(non_roots)  # number of non-root vertices
    e = len(edges)      # number of edges (= n - |roots| = m)
    
    B = matrix(ZZ, m, e)
    
    for j, (u, v) in enumerate(edges):
        # Orient edge from u to v (lower to higher)
        if u > v:
            u, v = v, u
        
        # u -> v means: -1 at u, +1 at v
        if u in non_roots:
            i = non_roots.index(u)
            B[i, j] = -1
        if v in non_roots:
            i = non_roots.index(v)
            B[i, j] = +1
    
    return B


def hypothesis_incidence_det_sign(forest, n, roots):
    """
    Hypothesis: sign = sign of determinant of incidence matrix.
    
    The incidence matrix is square (3x3 for n=6, roots={0,1,2}).
    """
    B = compute_incidence_matrix(forest, n, roots)
    
    if B.nrows() != B.ncols():
        return None
    
    det_B = B.det()
    
    if det_B > 0:
        return 1
    elif det_B < 0:
        return -1
    else:
        return 0


# ============================================================
# HYPOTHESIS 5: Product over edge pairs
# ============================================================

def hypothesis_edge_pair_product(forest, n):
    """
    Hypothesis: sign involves parity of some edge pair structure.
    
    Consider all pairs of edges and check some combinatorial property.
    """
    edges = list(forest)
    
    # Count pairs where edges "cross" in some sense
    crossings = 0
    for i in range(len(edges)):
        for j in range(i+1, len(edges)):
            e1 = edges[i]
            e2 = edges[j]
            # Crossing: min(e1) < min(e2) < max(e1) < max(e2) or vice versa
            a, b = min(e1), max(e1)
            c, d = min(e2), max(e2)
            if (a < c < b < d) or (c < a < d < b):
                crossings += 1
    
    return (-1)**crossings


# ============================================================
# MAIN TEST: Find which hypothesis matches
# ============================================================

def test_hypotheses(num_samples=10):
    """
    Test all hypotheses against actual forest signs.
    """
    print("=" * 70)
    print("TESTING FOREST SIGN HYPOTHESES")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    
    # Enumerate all 108 forests once
    all_forests = list(enumerate_rooted_forests(n, roots))
    print(f"Total forests: {len(all_forests)}")
    
    # Track match rates for each hypothesis
    hypotheses = {
        'crossing_012_345': lambda f: hypothesis_crossing_parity(f, ((0,1,2), (3,4,5))),
        'crossing_01_2345': lambda f: hypothesis_crossing_parity(f, ((0,1), (2,3,4,5))),
        'edge_orientation': hypothesis_edge_orientation,
        'inversions': lambda f: hypothesis_inversion_parity(f, n, roots),
        'incidence_det': lambda f: hypothesis_incidence_det_sign(f, n, roots),
        'edge_crossings': lambda f: hypothesis_edge_pair_product(f, n),
    }
    
    match_counts = {name: [] for name in hypotheses}
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 37)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        # Compute actual signs for all forests
        actual_signs = {}
        for forest in all_forests:
            sign, weight = compute_forest_sign(forest, lambdas, tilde_lambdas, x_spinor, y_spinor)
            if sign is not None:
                actual_signs[forest] = sign
        
        # Test each hypothesis
        for name, hyp_func in hypotheses.items():
            matches = 0
            total = 0
            for forest in all_forests:
                if forest not in actual_signs:
                    continue
                
                predicted = hyp_func(forest)
                actual = actual_signs[forest]
                
                if predicted == actual:
                    matches += 1
                total += 1
            
            if total > 0:
                match_counts[name].append(float(matches) / total)
    
    # Report results
    print("\nHypothesis Match Rates:")
    print("-" * 50)
    
    for name, rates in sorted(match_counts.items(), key=lambda x: -sum(x[1])/len(x[1]) if x[1] else 0):
        if rates:
            avg_rate = sum(rates) / len(rates)
            min_rate = min(rates)
            max_rate = max(rates)
            
            # A perfect match would be 100% every time
            # A anti-correlation would be 0% (flip and it's 100%)
            # Random would be ~50%
            
            status = ""
            if avg_rate > 0.99:
                status = "✓ PERFECT MATCH"
            elif avg_rate < 0.01:
                status = "✓ PERFECT ANTI-MATCH (flip sign)"
            elif avg_rate > 0.9:
                status = "~ Strong match"
            elif avg_rate < 0.1:
                status = "~ Strong anti-match"
            else:
                status = "✗ No correlation"
            
            print(f"  {name}: {avg_rate:.1%} (range: {min_rate:.1%}-{max_rate:.1%}) {status}")


def analyze_incidence_determinant():
    """
    Deep dive into the incidence matrix determinant hypothesis.
    
    The incidence matrix of a forest with 3 edges and 3 non-root vertices
    is a 3x3 matrix. Its determinant is ±1 for spanning forests.
    """
    print("\n" + "=" * 70)
    print("DETAILED ANALYSIS: Incidence Matrix Determinant")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    
    det_counts = {}
    
    for forest in all_forests:
        B = compute_incidence_matrix(forest, n, roots)
        det_B = B.det()
        det_counts[det_B] = det_counts.get(det_B, 0) + 1
    
    print(f"\nDeterminant distribution ({len(all_forests)} forests):")
    for det_val, count in sorted(det_counts.items()):
        pct = 100.0 * count / len(all_forests)
        print(f"  det = {det_val:+d}: {count} forests ({pct:.1f}%)")
    
    # If determinant is always ±1, then det sign is a candidate
    if set(det_counts.keys()) == {-1, 1}:
        print("\n✓ Determinant is always ±1!")
        print(f"  +1: {det_counts.get(1, 0)} forests")
        print(f"  -1: {det_counts.get(-1, 0)} forests")
        
        if det_counts.get(1, 0) == det_counts.get(-1, 0):
            print("  → Perfect 50/50 split!")
        
        return True
    
    return False


def discover_combinatorial_rule():
    """
    Try to discover a purely combinatorial rule for the sign.
    
    Strategy: Look for patterns that don't depend on kinematics.
    """
    print("\n" + "=" * 70)
    print("COMBINATORIAL SIGN RULE DISCOVERY")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    non_roots = [3, 4, 5]
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    
    # For each forest, compute various invariants
    forest_data = []
    
    for forest in all_forests:
        # Basic properties
        edges = list(forest)
        
        # Incidence matrix determinant
        B = compute_incidence_matrix(forest, n, roots)
        det_B = B.det()
        
        # Edge crossings
        crossings = 0
        for i in range(len(edges)):
            for j in range(i+1, len(edges)):
                e1, e2 = edges[i], edges[j]
                a, b = min(e1), max(e1)
                c, d = min(e2), max(e2)
                if (a < c < b < d) or (c < a < d < b):
                    crossings += 1
        
        # Edges involving specific vertices
        edges_with_3 = sum(1 for e in edges if 3 in e)
        edges_with_4 = sum(1 for e in edges if 4 in e)
        edges_with_5 = sum(1 for e in edges if 5 in e)
        
        # Edges between non-roots
        edges_between_non_roots = sum(1 for e in edges 
                                       if e[0] in non_roots and e[1] in non_roots)
        
        forest_data.append({
            'forest': forest,
            'det_B': det_B,
            'crossings': crossings,
            'edges_with_3': edges_with_3,
            'edges_with_4': edges_with_4,
            'edges_with_5': edges_with_5,
            'edges_between_non_roots': edges_between_non_roots,
        })
    
    # Now test if kinematics-dependent sign matches any combinatorial invariant
    print("\nTesting incidence determinant as the sign rule...")
    
    result = sample_spinor_helicity_conserving(n=6, seed=42)
    if result is None:
        print("Sampling failed")
        return
    
    lambdas, tilde_lambdas = result
    x_spinor = vector(QQ, [1, 2])
    y_spinor = vector(QQ, [3, 1])
    
    matches = 0
    anti_matches = 0
    other = 0
    
    for fd in forest_data:
        forest = fd['forest']
        sign_actual, _ = compute_forest_sign(forest, lambdas, tilde_lambdas, x_spinor, y_spinor)
        
        if sign_actual is None:
            continue
        
        sign_det = fd['det_B']
        
        if sign_actual == sign_det:
            matches += 1
        elif sign_actual == -sign_det:
            anti_matches += 1
        else:
            other += 1
    
    total = matches + anti_matches + other
    print(f"\nIncidence det vs actual sign:")
    print(f"  Matches: {matches}/{total} = {100.0*matches/total:.1f}%")
    print(f"  Anti-matches: {anti_matches}/{total} = {100.0*anti_matches/total:.1f}%")
    
    # Check multiple samples
    print("\nTesting across multiple kinematic samples...")
    
    det_match_rates = []
    
    for seed in range(20):
        result = sample_spinor_helicity_conserving(n=6, seed=seed * 17)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        matches = 0
        total = 0
        
        for fd in forest_data:
            forest = fd['forest']
            sign_actual, _ = compute_forest_sign(forest, lambdas, tilde_lambdas, x_spinor, y_spinor)
            
            if sign_actual is None:
                continue
            
            total += 1
            if sign_actual == fd['det_B']:
                matches += 1
        
        if total > 0:
            det_match_rates.append(float(matches) / total)
    
    if det_match_rates:
        avg = sum(det_match_rates) / len(det_match_rates)
        print(f"\nAverage det match rate: {avg:.1%}")
        
        if avg > 0.99:
            print("✓ INCIDENCE DETERMINANT IS THE SIGN RULE!")
        elif avg < 0.01:
            print("✓ NEGATIVE INCIDENCE DETERMINANT IS THE SIGN RULE!")
        else:
            print("✗ Sign depends on kinematics, not just combinatorics")


if __name__ == "__main__":
    # First, analyze the incidence determinant structure
    analyze_incidence_determinant()
    
    # Then test all hypotheses
    test_hypotheses(num_samples=10)
    
    # Deep dive into combinatorial rule
    discover_combinatorial_rule()

