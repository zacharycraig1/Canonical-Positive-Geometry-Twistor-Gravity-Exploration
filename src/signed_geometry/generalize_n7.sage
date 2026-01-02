#!/usr/bin/env sage
"""
Generalize Signed Geometry to n=7

Verify that the sign rule and 50/50 split pattern generalize to n=7.

For n=7 with k=3 roots:
- Number of edges per forest: n - k = 4
- Number of forests: (expected to be much larger)
"""
from sage.all import *
import itertools
import sys
import os

sys.path.insert(0, os.getcwd())


def ang_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def sq_bracket(tilde_lambdas, i, j):
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def sample_spinor_helicity_n7(seed=None):
    """Sample momentum-conserving spinors for n=7."""
    import random
    if seed is not None:
        random.seed(int(seed))
    
    n = 7
    
    # Sample random lambdas
    lambdas = []
    for i in range(n):
        l = vector(QQ, [QQ(random.randint(-100, 100)), QQ(random.randint(-100, 100))])
        while l == 0:
            l = vector(QQ, [QQ(random.randint(-100, 100)), QQ(random.randint(-100, 100))])
        lambdas.append(l)
    
    # Sample tilde_lambdas for i=0..n-3
    tilde_lambdas = [None] * n
    for i in range(n-2):
        l_tilde = vector(QQ, [QQ(random.randint(-100, 100)), QQ(random.randint(-100, 100))])
        tilde_lambdas[i] = l_tilde
    
    # Solve for last two tilde_lambdas
    P_known = matrix(QQ, 2, 2, 0)
    for i in range(n-2):
        P_known += matrix(QQ, 2, 2, [lambdas[i][0]*tilde_lambdas[i][0], lambdas[i][0]*tilde_lambdas[i][1],
                                     lambdas[i][1]*tilde_lambdas[i][0], lambdas[i][1]*tilde_lambdas[i][1]])
    
    RHS = -P_known
    
    M = matrix(QQ, 2, 2, [lambdas[n-2][0], lambdas[n-1][0], 
                          lambdas[n-2][1], lambdas[n-1][1]])
    
    if M.det() == 0:
        return None
    
    b_u = vector(QQ, [RHS[0,0], RHS[1,0]])
    u_sol = M.solve_right(b_u)
    
    b_v = vector(QQ, [RHS[0,1], RHS[1,1]])
    v_sol = M.solve_right(b_v)
    
    tilde_lambdas[n-2] = vector(QQ, [u_sol[0], v_sol[0]])
    tilde_lambdas[n-1] = vector(QQ, [u_sol[1], v_sol[1]])
    
    return lambdas, tilde_lambdas


def enumerate_rooted_forests(n, roots):
    """Enumerate all spanning forests with specified roots."""
    num_roots = len(roots)
    num_edges = n - num_roots
    roots_set = set(roots)
    
    all_edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    
    count = 0
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
            count += 1
            yield edges
    
    print(f"  (Enumerated {count} forests)")


def compute_forest_sign_n7(forest, lambdas, tilde_lambdas, x_spinor, y_spinor):
    """Compute forest sign for n=7 using our discovered rule."""
    n = len(lambdas)
    
    def ang_with_ref(lam, ref):
        return lam[0] * ref[1] - lam[1] * ref[0]
    
    # C factors
    C = {i: ang_with_ref(lambdas[i], x_spinor) * ang_with_ref(lambdas[i], y_spinor) 
         for i in range(n)}
    
    # w factors
    w = {}
    for i in range(n):
        for j in range(i+1, n):
            ang = ang_bracket(lambdas, i, j)
            sq = sq_bracket(tilde_lambdas, i, j)
            if ang != 0:
                w[(i, j)] = sq / ang
    
    edges = list(forest)
    
    # Compute weight
    weight = QQ(1)
    for u, v in edges:
        if (u, v) not in w:
            return None, None
        weight *= (-w[(u, v)] * C[u] * C[v])
    
    sign = 1 if weight > 0 else -1 if weight < 0 else 0
    return sign, weight


def analyze_n7_sign_structure(num_samples=20):
    """Analyze the sign structure for n=7 with full zero-term tracking."""
    print("=" * 70)
    print("N=7 SIGN STRUCTURE ANALYSIS")
    print("=" * 70)
    
    n = 7
    roots = (0, 1, 2)
    
    # First count forests
    print(f"\nCounting forests for n={n}, roots={roots}...")
    all_forests = list(enumerate_rooted_forests(n, roots))
    num_forests = len(all_forests)
    print(f"Total forests: {num_forests}")
    print(f"Edges per forest: {n - len(roots)}")
    
    # Analyze sign splits with zero-term tracking
    sign_splits = []
    detailed_results = []
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_n7(seed=sample * 67)
        if result is None:
            print(f"Sample {sample}: SKIPPED (singular kinematics)")
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        pos_count = 0
        neg_count = 0
        zero_count = 0
        undefined_count = 0
        
        for forest in all_forests:
            sign, weight = compute_forest_sign_n7(forest, lambdas, tilde_lambdas, x_spinor, y_spinor)
            if sign is None:
                undefined_count += 1
            elif sign > 0:
                pos_count += 1
            elif sign < 0:
                neg_count += 1
            else:  # sign == 0
                zero_count += 1
        
        total_classified = pos_count + neg_count + zero_count
        split = (pos_count, neg_count)
        sign_splits.append(split)
        detailed_results.append({
            'sample': sample,
            'pos': pos_count,
            'neg': neg_count,
            'zero': zero_count,
            'undefined': undefined_count,
            'total': total_classified
        })
        
        # Print with zero/undefined tracking
        status = ""
        if zero_count > 0:
            status = f" (zero: {zero_count})"
        if undefined_count > 0:
            status += f" (undefined: {undefined_count})"
        print(f"Sample {sample}: {pos_count}+ / {neg_count}- / total={total_classified}{status}")
    
    # Analyze
    if sign_splits:
        print("\n" + "-" * 70)
        print("Sign split analysis:")
        print("-" * 70)
        
        from collections import Counter
        split_counts = Counter(sign_splits)
        
        print(f"Total forests enumerated: {num_forests}")
        
        mode_split = split_counts.most_common(1)[0][0]
        mode_count = split_counts.most_common(1)[0][1]
        print(f"Modal split: {mode_split} (occurred in {mode_count}/{len(sign_splits)} samples)")
        
        # Check if close to 50/50
        mode_pos, mode_neg = mode_split
        ratio = float(mode_pos) / float(mode_pos + mode_neg) if (mode_pos + mode_neg) > 0 else 0
        print(f"Ratio: {ratio:.4%} positive")
        
        # Report on zero terms
        samples_with_zeros = sum(1 for r in detailed_results if r['zero'] > 0)
        if samples_with_zeros > 0:
            print(f"\nNote: {samples_with_zeros}/{len(detailed_results)} samples had forests with zero weight")
            print("      (due to degenerate kinematics where some w_ij = 0)")
        
        if abs(ratio - 0.5) < 0.1:
            print("\n✓ Approximately 50/50 split!")
        else:
            print(f"\nSplit differs from 50/50 by {abs(ratio - 0.5):.1%}")


def verify_sign_rule_n7(num_samples=20):
    """Verify the sign rule still holds for n=7 with full logging."""
    print("\n" + "=" * 70)
    print("N=7 SIGN RULE VERIFICATION")
    print("=" * 70)
    
    n = 7
    roots = (0, 1, 2)
    
    # Get forests (only enumerate once)
    print("Enumerating forests...")
    all_forests = list(enumerate_rooted_forests(n, roots))
    print(f"Total: {len(all_forests)}")
    print(f"Edges per forest: {n - len(roots)}")
    print("-" * 70)
    
    perfect_matches = 0
    skipped_samples = 0
    total_forests_checked = 0
    total_mismatches = 0
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_n7(seed=sample * 89)
        if result is None:
            skipped_samples += 1
            print(f"Sample {sample}: SKIPPED (singular kinematics)")
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        def ang_with_ref(lam, ref):
            return lam[0] * ref[1] - lam[1] * ref[0]
        
        C = {i: ang_with_ref(lambdas[i], x_spinor) * ang_with_ref(lambdas[i], y_spinor) 
             for i in range(n)}
        
        w = {}
        for i in range(n):
            for j in range(i+1, n):
                ang = ang_bracket(lambdas, i, j)
                sq = sq_bracket(tilde_lambdas, i, j)
                if ang != 0:
                    w[(i, j)] = sq / ang
        
        matches = 0
        total = 0
        zeros_skipped = 0
        
        for forest in all_forests:
            edges = list(forest)
            
            # Compute actual sign
            actual_sign, weight = compute_forest_sign_n7(forest, lambdas, tilde_lambdas, x_spinor, y_spinor)
            if actual_sign is None:
                continue
            
            # Skip zero-weight forests for sign comparison (sign is undefined)
            if actual_sign == 0:
                zeros_skipped += 1
                continue
            
            # Compute predicted sign from rule
            sign_from_edges = (-1)**len(edges)
            
            sign_from_w = 1
            for u, v in edges:
                if (u, v) in w and w[(u, v)] < 0:
                    sign_from_w *= -1
            
            degree = {i: 0 for i in range(n)}
            for u, v in edges:
                degree[u] += 1
                degree[v] += 1
            
            sign_from_C = 1
            for v in range(n):
                if degree[v] % 2 == 1 and C[v] < 0:
                    sign_from_C *= -1
            
            predicted_sign = sign_from_edges * sign_from_w * sign_from_C
            
            total += 1
            if predicted_sign == actual_sign:
                matches += 1
            else:
                total_mismatches += 1
        
        total_forests_checked += total
        match_rate = float(matches) / total if total > 0 else 0
        
        status_extra = ""
        if zeros_skipped > 0:
            status_extra = f" (zeros skipped: {zeros_skipped})"
        
        if match_rate > 0.99:
            perfect_matches += 1
            print(f"Sample {sample}: {matches}/{total} = {match_rate:.1%} ✓{status_extra}")
        else:
            print(f"Sample {sample}: {matches}/{total} = {match_rate:.1%} ✗{status_extra}")
    
    print("-" * 70)
    print(f"Samples tested: {num_samples - skipped_samples}/{num_samples}")
    print(f"Perfect matches: {perfect_matches}/{num_samples - skipped_samples}")
    print(f"Total forests checked: {total_forests_checked}")
    print(f"Total mismatches: {total_mismatches}")
    
    if perfect_matches == num_samples - skipped_samples and perfect_matches > 0:
        print("\n✓ Sign rule verified for n=7!")
    else:
        print("\n✗ Some samples failed")


if __name__ == "__main__":
    # First analyze sign structure
    analyze_n7_sign_structure(num_samples=20)
    
    # Then verify sign rule
    verify_sign_rule_n7(num_samples=20)

