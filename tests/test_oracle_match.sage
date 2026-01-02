#!/usr/bin/env sage
"""
MTT Identity: Independent Implementation Check

This test verifies that the forest expansion matches the Laplacian determinant.
It provides an independent implementation check by comparing two different algorithms:

1. Oracle: Computes det(L^minor) directly via Sage's matrix determinant
2. Forest: Computes Σ_F Π_e (edge weights) via explicit enumeration

These two algorithms should give the same result by the Matrix-Tree Theorem.
This validates that our forest enumeration and weight computation are correctly
implemented.

Run with: sage tests/test_oracle_match.sage
"""
from sage.all import *
import itertools
import sys
import os

sys.path.insert(0, os.getcwd())


def ang_bracket(lambdas, i, j):
    """Angle bracket <ij>."""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def sq_bracket(tilde_lambdas, i, j):
    """Square bracket [ij]."""
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def build_weighted_laplacian(n, lambdas, tilde_lambdas, x_spinor, y_spinor):
    """
    Build the Hodges weighted Laplacian matrix.
    
    Returns L, C_values, w_values or (None, None, None) if singular.
    """
    def ang_with_ref(lam, ref):
        return lam[0] * ref[1] - lam[1] * ref[0]
    
    C = {}
    for i in range(n):
        C[i] = ang_with_ref(lambdas[i], x_spinor) * ang_with_ref(lambdas[i], y_spinor)
    
    w = {}
    for i in range(n):
        for j in range(i+1, n):
            ang = ang_bracket(lambdas, i, j)
            sq = sq_bracket(tilde_lambdas, i, j)
            if ang == 0:
                return None, None, None
            w[(i, j)] = sq / ang
            w[(j, i)] = w[(i, j)]
    
    L = matrix(QQ, n, n)
    for i in range(n):
        row_sum = QQ(0)
        for j in range(n):
            if i != j:
                L[i, j] = -w[(i, j)] * C[i] * C[j]
                row_sum += w[(i, j)] * C[i] * C[j]
        L[i, i] = row_sum
    
    return L, C, w


def enumerate_rooted_forests(n, roots):
    """Enumerate all spanning forests of K_n with roots R."""
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


def compute_forest_sum_mtt(forests, w, C):
    """
    Compute forest sum using STANDARD MTT weights: a_ij = w_ij * C_i * C_j
    
    This should equal det(L^minor) exactly.
    """
    total = QQ(0)
    for forest in forests:
        term = QQ(1)
        for u, v in forest:
            if u > v:
                u, v = v, u
            term *= w[(u, v)] * C[u] * C[v]
        total += term
    return total


def compute_forest_sum_signed(forests, w, C):
    """
    Compute forest sum using SIGNED weights: b_ij = -w_ij * C_i * C_j
    
    This equals (-1)^|E| * det(L^minor).
    """
    total = QQ(0)
    for forest in forests:
        term = QQ(1)
        for u, v in forest:
            if u > v:
                u, v = v, u
            term *= (-w[(u, v)] * C[u] * C[v])
        total += term
    return total


def test_mtt_identity(num_samples=50):
    """
    TEST 1: Matrix-Tree Theorem Identity
    
    Verify: det(L^minor) = Σ_F Π_e a_e  (using standard MTT weights)
    
    This is the fundamental identity that our sign rule depends on.
    """
    print("=" * 70)
    print("TEST 1: Matrix-Tree Theorem Identity (non-circular)")
    print("=" * 70)
    print("Oracle: det(L^minor) via Sage matrix determinant")
    print("Forest: Σ_F Π_e (w_e × C_i × C_j) via explicit enumeration")
    print("-" * 70)
    
    load('src/spinor_sampling.sage')
    
    n = 6
    roots = (0, 1, 2)
    non_roots = [i for i in range(n) if i not in roots]
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    print(f"Total forests: {len(all_forests)}")
    
    matches = 0
    skipped = 0
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 37 + 11)
        if result is None:
            skipped += 1
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        L, C, w = build_weighted_laplacian(n, lambdas, tilde_lambdas, x_spinor, y_spinor)
        if L is None:
            skipped += 1
            continue
        
        # ORACLE: Direct determinant computation
        L_minor = L.matrix_from_rows_and_columns(non_roots, non_roots)
        det_oracle = L_minor.det()
        
        # FOREST: Explicit sum using standard MTT weights
        det_forest = compute_forest_sum_mtt(all_forests, w, C)
        
        # Compare
        if det_oracle == det_forest:
            matches += 1
        else:
            print(f"  MISMATCH at sample {sample}!")
            print(f"    Oracle: {det_oracle}")
            print(f"    Forest: {det_forest}")
            print(f"    Ratio:  {det_oracle / det_forest if det_forest != 0 else 'N/A'}")
    
    tested = num_samples - skipped
    print(f"\nResults: {matches}/{tested} exact matches ({skipped} skipped)")
    
    if matches == tested and tested > 0:
        print("✓ PASS: Matrix-Tree Theorem identity verified!")
        return True
    else:
        print("✗ FAIL: Matrix-Tree Theorem identity failed!")
        return False


def test_sign_convention_relation(num_samples=50):
    """
    TEST 2: Sign Convention Relation
    
    Verify: Σ_F Π_e b_e = (-1)^|E| × det(L^minor)
    
    where b_e = -a_e (signed-edge weights).
    """
    print("\n" + "=" * 70)
    print("TEST 2: Sign Convention Relation")
    print("=" * 70)
    print("Verify: signed_sum = (-1)^|E| × det(L^minor)")
    print("-" * 70)
    
    load('src/spinor_sampling.sage')
    
    n = 6
    roots = (0, 1, 2)
    non_roots = [i for i in range(n) if i not in roots]
    num_edges = n - len(roots)  # = 3
    
    all_forests = list(enumerate_rooted_forests(n, roots))
    
    matches = 0
    skipped = 0
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 41 + 7)
        if result is None:
            skipped += 1
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        L, C, w = build_weighted_laplacian(n, lambdas, tilde_lambdas, x_spinor, y_spinor)
        if L is None:
            skipped += 1
            continue
        
        # Oracle
        L_minor = L.matrix_from_rows_and_columns(non_roots, non_roots)
        det_oracle = L_minor.det()
        
        # Signed forest sum
        signed_sum = compute_forest_sum_signed(all_forests, w, C)
        
        # Expected relation
        expected = (-1)**num_edges * det_oracle
        
        if signed_sum == expected:
            matches += 1
        else:
            print(f"  MISMATCH at sample {sample}!")
    
    tested = num_samples - skipped
    print(f"\nResults: {matches}/{tested} exact matches ({skipped} skipped)")
    
    if matches == tested and tested > 0:
        print("✓ PASS: Sign convention relation verified!")
        return True
    else:
        print("✗ FAIL: Sign convention relation failed!")
        return False


def test_reference_independence(num_samples=20):
    """
    TEST 3: Reference Spinor Independence
    
    The full amplitude should be independent of the choice of reference spinors (x, y).
    """
    print("\n" + "=" * 70)
    print("TEST 3: Reference Spinor Independence")
    print("=" * 70)
    print("Full amplitude should not depend on choice of (x, y)")
    print("-" * 70)
    
    load('src/spinor_sampling.sage')
    
    n = 6
    roots = (0, 1, 2)
    non_roots = [i for i in range(n) if i not in roots]
    
    # Different reference spinor pairs
    ref_pairs = [
        (vector(QQ, [1, 0]), vector(QQ, [0, 1])),
        (vector(QQ, [1, 2]), vector(QQ, [3, 1])),
        (vector(QQ, [2, 5]), vector(QQ, [7, 3])),
        (vector(QQ, [1, 1]), vector(QQ, [1, -1])),
    ]
    
    consistent = 0
    skipped = 0
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 53 + 13)
        if result is None:
            skipped += 1
            continue
        
        lambdas, tilde_lambdas = result
        
        amplitudes = []
        valid = True
        
        for x_spinor, y_spinor in ref_pairs:
            L, C, w = build_weighted_laplacian(n, lambdas, tilde_lambdas, x_spinor, y_spinor)
            if L is None:
                valid = False
                break
            
            # Compute full amplitude with normalization
            L_minor = L.matrix_from_rows_and_columns(non_roots, non_roots)
            det_L = L_minor.det()
            
            # Helicity factor
            ang_01 = ang_bracket(lambdas, 0, 1)
            h_factor = ang_01**8
            
            # Normalization
            ang_12 = ang_bracket(lambdas, roots[0], roots[1])
            ang_23 = ang_bracket(lambdas, roots[1], roots[2])
            ang_31 = ang_bracket(lambdas, roots[2], roots[0])
            norm_R = (ang_12 * ang_23 * ang_31)**2
            
            # C factors for non-roots
            def ang_with_ref(lam, ref):
                return lam[0] * ref[1] - lam[1] * ref[0]
            
            prod_C_sq = QQ(1)
            for k in non_roots:
                C_k = ang_with_ref(lambdas[k], x_spinor) * ang_with_ref(lambdas[k], y_spinor)
                if C_k == 0:
                    valid = False
                    break
                prod_C_sq *= C_k**2
            
            if not valid or norm_R == 0:
                valid = False
                break
            
            # Full amplitude
            sign_factor = (-1)**(n - 1)
            amp = sign_factor * h_factor * det_L / (norm_R * prod_C_sq)
            amplitudes.append(amp)
        
        if not valid or len(amplitudes) < len(ref_pairs):
            skipped += 1
            continue
        
        # Check all amplitudes are equal
        if all(a == amplitudes[0] for a in amplitudes):
            consistent += 1
        else:
            print(f"  Sample {sample}: Reference dependence detected!")
            for i, (ref, amp) in enumerate(zip(ref_pairs, amplitudes)):
                print(f"    Ref {i}: {float(amp):.6e}")
    
    tested = num_samples - skipped
    print(f"\nResults: {consistent}/{tested} consistent across references ({skipped} skipped)")
    
    if consistent == tested and tested > 0:
        print("✓ PASS: Reference spinor independence verified!")
        return True
    else:
        print("✗ FAIL: Reference spinor dependence detected!")
        return False


def run_all_tests():
    """Run all MTT consistency tests."""
    print("=" * 70)
    print("MTT ↔ det(L) CONSISTENCY CHECK SUITE")
    print("=" * 70)
    print("These tests verify the forest expansion against Laplacian determinant.")
    print("=" * 70)
    
    results = {}
    
    results['mtt_identity'] = test_mtt_identity(num_samples=50)
    results['sign_convention'] = test_sign_convention_relation(num_samples=50)
    results['reference_independence'] = test_reference_independence(num_samples=20)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    all_pass = True
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_pass = False
    
    print("-" * 70)
    if all_pass:
        print("✓ ALL MTT CONSISTENCY TESTS PASSED")
    else:
        print("✗ SOME TESTS FAILED")
    
    return all_pass


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

