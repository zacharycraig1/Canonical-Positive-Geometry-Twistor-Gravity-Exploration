#!/usr/bin/env sage
"""
Kinematic Sign Analysis for Forest Terms

Since the sign is NOT purely combinatorial, it must depend on the
kinematic weights w_ij = [ij]/⟨ij⟩.

Key insight: Each forest term has the form:
    weight(F) = prod_{(i,j) in E(F)} (-w_ij * C_i * C_j)

The sign depends on:
1. The number of edges (always 3 for n=6, so (-1)^3 = -1 factor)
2. The signs of individual w_ij
3. The signs of individual C_i factors (squared, so always positive)

Wait - C_i appears in multiple edges. Need to track carefully.
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
    """Enumerate all spanning forests with specified roots."""
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


def analyze_forest_sign_structure(lambdas, tilde_lambdas, x_spinor, y_spinor, roots=(0, 1, 2)):
    """
    Decompose the sign of each forest term into its components.
    
    For a forest F with edges E(F):
        weight(F) = prod_{(i,j) in E(F)} (-w_ij * C_i * C_j)
        
    This factors as:
        weight(F) = (-1)^|E(F)| * prod_e w_e * prod_v C_v^{deg(v)}
        
    For n=6, k=3 roots:
        |E(F)| = 3 (always), so (-1)^3 = -1
        
    The sign is:
        sign(F) = -1 * sign(prod_e w_e) * sign(prod_v C_v^{deg(v)})
    """
    n = len(lambdas)
    
    def ang_with_ref(lam, ref):
        return lam[0] * ref[1] - lam[1] * ref[0]
    
    # Compute C_i values
    C = {}
    for i in range(n):
        C[i] = ang_with_ref(lambdas[i], x_spinor) * ang_with_ref(lambdas[i], y_spinor)
    
    # Compute w_ij values
    w = {}
    for i in range(n):
        for j in range(i+1, n):
            ang = ang_bracket(lambdas, i, j)
            sq = sq_bracket(tilde_lambdas, i, j)
            if ang != 0:
                w[(i, j)] = sq / ang
            else:
                w[(i, j)] = None
    
    # Analyze each forest
    forests = list(enumerate_rooted_forests(n, roots))
    forest_analysis = []
    
    for forest in forests:
        edges = list(forest)
        
        # Count vertex degrees
        degree = {i: 0 for i in range(n)}
        for u, v in edges:
            degree[u] += 1
            degree[v] += 1
        
        # Compute product of w_ij
        w_product = QQ(1)
        w_signs = []
        valid = True
        for u, v in edges:
            if w[(u, v)] is None:
                valid = False
                break
            w_product *= w[(u, v)]
            w_signs.append(1 if w[(u, v)] > 0 else -1)
        
        if not valid:
            continue
        
        # Compute product of C_v^{deg(v)}
        C_product = QQ(1)
        C_signs = []
        for v in range(n):
            if degree[v] > 0:
                C_product *= C[v] ** degree[v]
                C_signs.append((v, degree[v], 1 if C[v] > 0 else -1))
        
        # Total sign
        total_weight = (-1)**len(edges) * w_product * C_product
        sign = 1 if total_weight > 0 else -1
        
        # Sign decomposition
        sign_from_minus_one = (-1)**len(edges)
        sign_from_w = 1 if w_product > 0 else -1
        sign_from_C = 1 if C_product > 0 else -1
        
        forest_analysis.append({
            'forest': forest,
            'sign': sign,
            'sign_from_(-1)^|E|': sign_from_minus_one,
            'sign_from_w': sign_from_w,
            'sign_from_C': sign_from_C,
            'w_signs': tuple(w_signs),
            'degrees': tuple(degree[i] for i in range(n)),
        })
    
    return forest_analysis, C, w


def analyze_sign_factorization(num_samples=5):
    """
    Analyze how the sign factors across multiple kinematic samples.
    """
    print("=" * 70)
    print("SIGN FACTORIZATION ANALYSIS")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    print(f"Total forests: {len(all_forests)}")
    
    # Track which components determine the sign
    for sample in range(num_samples):
        print(f"\n--- Sample {sample} ---")
        
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 53)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        analysis, C, w = analyze_forest_sign_structure(
            lambdas, tilde_lambdas, x_spinor, y_spinor, roots
        )
        
        # Print C signs
        C_signs = {i: (1 if C[i] > 0 else -1) for i in range(n)}
        print(f"C signs: {C_signs}")
        
        # Print w signs
        w_pos = sum(1 for (i,j), val in w.items() if val is not None and val > 0)
        w_neg = sum(1 for (i,j), val in w.items() if val is not None and val < 0)
        print(f"w signs: {w_pos} positive, {w_neg} negative")
        
        # Sign distribution
        pos_forests = sum(1 for fa in analysis if fa['sign'] > 0)
        neg_forests = sum(1 for fa in analysis if fa['sign'] < 0)
        print(f"Forest signs: {pos_forests}+, {neg_forests}-")
        
        # Which factor dominates?
        # Group forests by their w_sign pattern
        w_pattern_counts = {}
        for fa in analysis:
            pattern = fa['w_signs']
            if pattern not in w_pattern_counts:
                w_pattern_counts[pattern] = {'pos': 0, 'neg': 0}
            if fa['sign'] > 0:
                w_pattern_counts[pattern]['pos'] += 1
            else:
                w_pattern_counts[pattern]['neg'] += 1


def find_kinematic_sign_rule():
    """
    The sign depends on the signs of w_ij.
    
    Hypothesis: sign(F) = product of edge contributions, where each edge
    contributes based on the sign of w_ij.
    
    Since weight(F) = prod_e (-w_e * C^2_e), and C^2 > 0, the sign is:
    
        sign(F) = prod_e sign(-w_e) = (-1)^|E| * prod_e sign(w_e)
        
    For |E| = 3:
        sign(F) = -prod_e sign(w_e)
    
    Let's verify this!
    """
    print("\n" + "=" * 70)
    print("KINEMATIC SIGN RULE VERIFICATION")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    
    perfect_matches = 0
    total_tests = 0
    
    for sample in range(20):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 41)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        analysis, C, w = analyze_forest_sign_structure(
            lambdas, tilde_lambdas, x_spinor, y_spinor, roots
        )
        
        matches = 0
        for fa in analysis:
            # Predicted sign from w signs only
            # Each term is -w_ij * (C factors)^2
            # C factors are squared, so always positive
            # Actually wait - C_i * C_j, not C_i^2
            
            # Let's think again...
            # For edge (i,j): contribution is -w_ij * C_i * C_j
            # 
            # If we factor out C completely:
            # weight = prod_e (-w_e * C_{e_1} * C_{e_2})
            #        = (-1)^|E| * prod_e w_e * prod_e C_{e_1} * C_{e_2}
            #        = (-1)^|E| * prod_e w_e * prod_v C_v^{deg(v)}
            #
            # The sign is:
            # sign = (-1)^|E| * sign(prod_e w_e) * sign(prod_v C_v^{deg(v)})
            #
            # For a spanning forest: sum of degrees = 2|E| (each edge has 2 endpoints)
            # And product C_v^{deg(v)} = product over vertices of C_v raised to degree
            
            # But in the formula, it's C_i * C_j per edge, not squared.
            # So we have prod_e (C_{e_1} * C_{e_2}).
            # 
            # A vertex of degree d appears in d edges.
            # So the total power of C_v is d.
            # product = prod_v C_v^{deg(v)}
            
            forest = fa['forest']
            edges = list(forest)
            
            # Count degrees
            degree = {i: 0 for i in range(n)}
            for u, v in edges:
                degree[u] += 1
                degree[v] += 1
            
            # Sign from C factors
            sign_C = 1
            for v in range(n):
                if degree[v] % 2 == 1:  # Odd power of C_v
                    if C[v] < 0:
                        sign_C *= -1
            
            # Sign from w factors
            sign_w = 1
            for u, v in edges:
                if w[(u, v)] < 0:
                    sign_w *= -1
            
            # Total predicted sign
            sign_from_edges = (-1)**len(edges)  # = -1 for 3 edges
            predicted_sign = sign_from_edges * sign_w * sign_C
            
            actual_sign = fa['sign']
            
            if predicted_sign == actual_sign:
                matches += 1
        
        match_rate = float(matches) / len(analysis) if analysis else 0
        
        if match_rate > 0.99:
            perfect_matches += 1
        
        total_tests += 1
    
    print(f"\nPerfect match rate: {perfect_matches}/{total_tests}")
    
    if perfect_matches == total_tests:
        print("\n✓ SIGN RULE DISCOVERED!")
        print("\nThe sign of forest F is:")
        print("  sign(F) = (-1)^|E| × sign(∏_e w_e) × sign(∏_v C_v^{deg(v)})")
        print("\nWhere:")
        print("  - |E| = number of edges (3 for n=6)")
        print("  - w_e = [ij]/⟨ij⟩ for edge (i,j)")
        print("  - C_v = ⟨v,x⟩⟨v,y⟩ for reference spinors x,y")
        print("  - deg(v) = degree of vertex v in the forest")


def analyze_C_dependence():
    """
    Since C_i depends on reference spinors, and the amplitude is
    reference-spinor-independent, the sign pattern must also be independent.
    
    Let's verify that the overall sign pattern is independent of reference choice.
    """
    print("\n" + "=" * 70)
    print("REFERENCE SPINOR INDEPENDENCE CHECK")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    
    result = sample_spinor_helicity_conserving(n=6, seed=12345)
    if result is None:
        print("Sampling failed")
        return
    
    lambdas, tilde_lambdas = result
    
    # Try different reference spinors
    reference_pairs = [
        (vector(QQ, [1, 2]), vector(QQ, [3, 1])),
        (vector(QQ, [1, 0]), vector(QQ, [0, 1])),
        (vector(QQ, [5, 7]), vector(QQ, [11, 3])),
        (vector(QQ, [-2, 5]), vector(QQ, [7, -3])),
    ]
    
    sign_patterns = []
    
    for x_spinor, y_spinor in reference_pairs:
        analysis, C, w = analyze_forest_sign_structure(
            lambdas, tilde_lambdas, x_spinor, y_spinor, roots
        )
        
        signs = tuple(fa['sign'] for fa in sorted(analysis, key=lambda x: x['forest']))
        sign_patterns.append(signs)
        
        pos = sum(1 for s in signs if s > 0)
        neg = sum(1 for s in signs if s < 0)
        print(f"Ref ({x_spinor}, {y_spinor}): {pos}+, {neg}-")
    
    # Check if patterns are the same or related
    if len(set(sign_patterns)) == 1:
        print("\n✓ Sign pattern is IDENTICAL for all reference spinors!")
    else:
        print("\n✗ Sign pattern varies with reference spinors")
        
        # Check if they're related by overall sign flip
        base = sign_patterns[0]
        for i, pattern in enumerate(sign_patterns[1:], 1):
            if all(s1 == s2 for s1, s2 in zip(base, pattern)):
                print(f"  Pattern {i+1} = Pattern 1")
            elif all(s1 == -s2 for s1, s2 in zip(base, pattern)):
                print(f"  Pattern {i+1} = -Pattern 1 (overall flip)")
            else:
                matching = sum(1 for s1, s2 in zip(base, pattern) if s1 == s2)
                print(f"  Pattern {i+1}: {matching}/{len(base)} match Pattern 1")


if __name__ == "__main__":
    analyze_sign_factorization(num_samples=3)
    find_kinematic_sign_rule()
    analyze_C_dependence()

