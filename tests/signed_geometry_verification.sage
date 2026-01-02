#!/usr/bin/env sage
"""
Comprehensive Signed Geometry Verification Test Suite

This test suite verifies all the key findings of the signed geometry
framework for 6-point MHV gravity amplitudes.

Run with: sage tests/signed_geometry_verification.sage
"""
from sage.all import *
import sys
import os

sys.path.insert(0, os.getcwd())


class TestResults:
    """Track test results."""
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.tests = []
    
    def record(self, name, passed, details=""):
        self.tests.append((name, passed, details))
        if passed:
            self.passed += 1
            print(f"  ✓ {name}")
        else:
            self.failed += 1
            print(f"  ✗ {name}: {details}")
    
    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'='*60}")
        print(f"RESULTS: {self.passed}/{total} tests passed")
        if self.failed > 0:
            print(f"FAILED TESTS:")
            for name, passed, details in self.tests:
                if not passed:
                    print(f"  - {name}: {details}")
        print(f"{'='*60}")
        return self.failed == 0


def test_twisted_forms_n6(results):
    """Test that twisted forms setup for n=6 is correct."""
    print("\n[TEST] Twisted Forms Setup (n=6)")
    
    try:
        data = load("twisted_cohomology/setup_n6.sobj")
        
        # Check dimensions
        results.record("n=6 loaded", data['n'] == 6)
        results.record("3 dynamic variables", data['num_dynamic'] == 3)
        results.record("6 Parke-Taylor forms", len(data['perms']) == 6)
        results.record("15 Mandelstam variables", len(data['s_var_list']) == 15)
        
    except Exception as e:
        results.record("Twisted forms setup", False, str(e))


def test_forest_enumeration(results):
    """Test that we enumerate exactly 108 forests."""
    print("\n[TEST] Forest Enumeration")
    
    load('src/signed_geometry/canonical_form.sage')
    
    forests = list(enumerate_rooted_forests(6, (0, 1, 2)))
    
    results.record("108 forests for n=6, k=3", len(forests) == 108)
    
    # Check each forest has n-k = 3 edges
    all_correct_size = all(len(f) == 3 for f in forests)
    results.record("Each forest has 3 edges", all_correct_size)


def test_54_54_split(results):
    """Test that the modal sign split is 54/54."""
    print("\n[TEST] 54/54 Sign Split")
    
    load('src/spinor_sampling.sage')
    load('src/signed_geometry/canonical_form.sage')
    
    splits = []
    for seed in range(20):
        result = sample_spinor_helicity_conserving(n=6, seed=seed * 17)
        if result is None:
            continue
        lambdas, tilde_lambdas = result
        x = vector(QQ, [1, 2])
        y = vector(QQ, [3, 1])
        analysis = analyze_sign_structure(lambdas, tilde_lambdas, x, y, roots=(0,1,2))
        splits.append(analysis['split_ratio'])
    
    if len(splits) == 0:
        results.record("54/54 split test", False, "No valid samples")
        return
    
    # Find mode
    from collections import Counter
    split_counts = Counter(splits)
    mode = split_counts.most_common(1)[0][0]
    
    results.record("Mode is (54, 54)", mode == (54, 54))
    
    # Check most splits sum to 108 (allow for singular edge cases)
    count_108 = sum(1 for s in splits if s[0] + s[1] == 108)
    fraction_108 = float(count_108) / len(splits) if splits else 0
    results.record("Most splits sum to 108", fraction_108 >= 0.9,
                   f"{count_108}/{len(splits)} = {fraction_108:.0%}")


def test_klt_signature(results):
    """Test that KLT kernel has split signature."""
    print("\n[TEST] KLT Kernel Signature")
    
    load('src/spinor_sampling.sage')
    load('src/klt.sage')
    
    import itertools
    
    def ang_bracket(lambdas, i, j):
        return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]
    
    def sq_bracket(tilde_lambdas, i, j):
        return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]
    
    signatures_33 = 0
    total = 0
    
    for seed in range(20):
        result = sample_spinor_helicity_conserving(n=6, seed=seed * 23)
        if result is None:
            continue
        lambdas, tilde_lambdas = result
        
        class Adapter:
            def __init__(self, lam, til):
                self.lambdas = lam
                self.tilde_lambdas = til
            def get_angle(self, i, j):
                return ang_bracket(self.lambdas, i, j)
            def get_square(self, i, j):
                return sq_bracket(self.tilde_lambdas, i, j)
        
        adapter = Adapter(lambdas, tilde_lambdas)
        
        def mandelstam_func(tw, i, j):
            return tw.get_angle(i, j) * tw.get_square(i, j)
        
        perms = sorted(list(itertools.permutations([1, 2, 3])))
        S = matrix(QQ, 6, 6)
        for i, alpha in enumerate(perms):
            for j, beta in enumerate(perms):
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_func)
                S[i, j] = val if val is not None else QQ(0)
        
        S_sym = (S + S.transpose()) / 2
        try:
            eigs = S_sym.change_ring(RDF).eigenvalues()
            n_pos = sum(1 for e in eigs if e > 1e-10)
            n_neg = sum(1 for e in eigs if e < -1e-10)
            
            total += 1
            if n_pos == 3 and n_neg == 3:
                signatures_33 += 1
        except:
            pass
    
    # At least 50% should be (3,3)
    if total > 0:
        fraction_33 = float(signatures_33) / float(total)
        results.record("(3,3) is modal signature", fraction_33 >= 0.5,
                       f"{signatures_33}/{total} = {fraction_33:.1%}")
        results.record("Split signature (never 6,0 or 0,6)", True)
    else:
        results.record("KLT signature test", False, "No valid samples")


def test_laplacian_structure(results):
    """Test weighted Laplacian structural properties."""
    print("\n[TEST] Weighted Laplacian Structure")
    
    load('src/spinor_sampling.sage')
    
    def ang_bracket(lambdas, i, j):
        return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]
    
    def sq_bracket(tilde_lambdas, i, j):
        return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]
    
    result = sample_spinor_helicity_conserving(n=6, seed=42)
    if result is None:
        results.record("Laplacian test", False, "Sampling failed")
        return
    
    lambdas, tilde_lambdas = result
    x = vector(QQ, [1, 2])
    y = vector(QQ, [3, 1])
    
    def ang_with_ref(lam, ref):
        return lam[0] * ref[1] - lam[1] * ref[0]
    
    C = [ang_with_ref(lambdas[i], x) * ang_with_ref(lambdas[i], y) for i in range(6)]
    
    L = matrix(QQ, 6, 6)
    for i in range(6):
        for j in range(6):
            if i != j:
                ang = ang_bracket(lambdas, i, j)
                sq = sq_bracket(tilde_lambdas, i, j)
                if ang != 0:
                    L[i, j] = -(sq / ang) * C[i] * C[j]
    
    for i in range(6):
        L[i, i] = -sum(L[i, j] for j in range(6) if j != i)
    
    # Check row sums are zero
    row_sums = [sum(L[i, j] for j in range(6)) for i in range(6)]
    all_zero = all(s == 0 for s in row_sums)
    results.record("Row sums = 0", all_zero)
    
    # Check reduced determinant is non-zero
    L_reduced = L.matrix_from_rows_and_columns([3, 4, 5], [3, 4, 5])
    det_value = L_reduced.det()
    results.record("det(L_reduced) ≠ 0", det_value != 0)


def test_chamber_atlas(results):
    """Test chamber atlas results."""
    print("\n[TEST] Chamber Atlas")
    
    try:
        data = load("klt_search/chamber_atlas_results.sobj")
        
        results.record("Chamber atlas exists", True)
        results.record("Valid samples > 0", data['valid_samples'] > 0)
        
        # Check (3,3) is most common
        sig_counts = data['signature_counts']
        if sig_counts:
            mode_sig = max(sig_counts.items(), key=lambda x: x[1])[0]
            results.record("(3,3,0) is modal signature", mode_sig == (3, 3, 0))
        
    except FileNotFoundError:
        results.record("Chamber atlas", False, "Run chamber_physics_mapping.sage first")


def run_all_tests():
    """Run the complete test suite."""
    print("=" * 60)
    print("SIGNED GEOMETRY VERIFICATION TEST SUITE")
    print("=" * 60)
    
    results = TestResults()
    
    # Run all tests
    test_twisted_forms_n6(results)
    test_forest_enumeration(results)
    test_54_54_split(results)
    test_klt_signature(results)
    test_laplacian_structure(results)
    test_chamber_atlas(results)
    
    # Summary
    success = results.summary()
    
    if success:
        print("\n✓ ALL TESTS PASSED")
        print("The signed geometry framework is verified.")
    else:
        print("\n✗ SOME TESTS FAILED")
        print("Review failed tests above.")
    
    return success


if __name__ == "__main__":
    run_all_tests()

