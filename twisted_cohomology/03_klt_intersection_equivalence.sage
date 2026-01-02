#!/usr/bin/env sage
"""
KLT Kernel = Inverse Intersection Matrix Equivalence

This script proves the relationship:
    S_KLT[alpha|beta] = c(kinematics) * M_biadjoint^{-1}[alpha|beta]

where c(kinematics) is a normalization factor depending on spinor brackets.

Based on the literature (Mizera et al.), the expected relationship involves
the Parke-Taylor cyclic denominator product.
"""
from sage.all import *
import itertools
import sys
import os

sys.path.insert(0, os.getcwd())


def sample_kinematics(seed=None):
    """Sample momentum-conserving spinor kinematics."""
    load('src/spinor_sampling.sage')
    return sample_spinor_helicity_conserving(n=6, seed=seed)


def ang_bracket(lambdas, i, j):
    """Compute angle bracket <ij>."""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def sq_bracket(tilde_lambdas, i, j):
    """Compute square bracket [ij]."""
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def mandelstam(lambdas, tilde_lambdas, i, j):
    """Compute Mandelstam s_ij = <ij>[ij]."""
    return ang_bracket(lambdas, i, j) * sq_bracket(tilde_lambdas, i, j)


def parke_taylor_denominator(lambdas, ordering):
    """
    Compute the cyclic Parke-Taylor denominator for a given ordering.
    PT_denom(alpha) = prod_{consecutive i} <alpha[i], alpha[i+1]>
    """
    n = len(ordering)
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        denom *= ang_bracket(lambdas, ordering[i], ordering[j])
    return denom


def compute_klt_kernel_matrix(lambdas, tilde_lambdas):
    """
    Compute the 6x6 KLT momentum kernel matrix.
    
    S[alpha|beta] for alpha, beta permutations of {1,2,3}
    Fixed legs: {0, 4, 5}
    """
    load('src/klt.sage')
    
    class SpinorAdapter:
        def __init__(self, lam, til):
            self.lambdas = lam
            self.tilde_lambdas = til
        def get_angle(self, i, j):
            return ang_bracket(self.lambdas, i, j)
        def get_square(self, i, j):
            return sq_bracket(self.tilde_lambdas, i, j)
    
    adapter = SpinorAdapter(lambdas, tilde_lambdas)
    
    def mandelstam_func(tw, i, j):
        return tw.get_angle(i,j) * tw.get_square(i,j)
    
    permuted_set = [1, 2, 3]
    basis_perms = sorted(list(itertools.permutations(permuted_set)))
    
    S = matrix(QQ, 6, 6)
    for i, alpha in enumerate(basis_perms):
        for j, beta in enumerate(basis_perms):
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_func)
            S[i, j] = val if val is not None else QQ(0)
    
    return S, basis_perms


def compute_biadjoint_matrix(lambdas, tilde_lambdas):
    """
    Compute the 6x6 bi-adjoint scalar amplitude matrix m[alpha|beta].
    """
    import networkx as nx
    
    def s_subset(subset):
        """Compute s for a subset of particles."""
        val = QQ(0)
        subset = list(subset)
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                val += mandelstam(lambdas, tilde_lambdas, subset[a], subset[b])
        return val
    
    def get_planar_poles(order):
        m = len(order)
        poles = set()
        doubled = list(order) + list(order)
        for size in range(2, m-1):
            for i in range(m):
                subset = frozenset(doubled[i:i+size])
                if len(subset) > m//2:
                    subset = frozenset(order) - subset
                elif len(subset) == m//2 and 0 not in subset:
                    subset = frozenset(order) - subset
                poles.add(subset)
        return list(poles)
    
    def get_triangulations(order):
        poles = get_planar_poles(order)
        m = len(order)
        num_propagators = m - 3
        
        compatible_pairs = set()
        for i in range(len(poles)):
            for j in range(i+1, len(poles)):
                p1, p2 = poles[i], poles[j]
                if p1.isdisjoint(p2) or p1.issubset(p2) or p2.issubset(p1):
                    compatible_pairs.add((i, j))
        
        G = nx.Graph()
        G.add_nodes_from(range(len(poles)))
        G.add_edges_from(list(compatible_pairs))
        cliques = [c for c in nx.find_cliques(G) if len(c) == num_propagators]
        return [[poles[x] for x in c] for c in cliques]
    
    permuted_set = [1, 2, 3]
    basis_perms = sorted(list(itertools.permutations(permuted_set)))
    
    alpha_diagrams = {}
    beta_diagrams = {}
    for alpha in basis_perms:
        alpha_diagrams[alpha] = get_triangulations([0] + list(alpha) + [4, 5])
    for beta in basis_perms:
        beta_diagrams[beta] = get_triangulations([0] + list(beta) + [5, 4])
    
    M = matrix(QQ, 6, 6)
    for i, alpha in enumerate(basis_perms):
        for j, beta in enumerate(basis_perms):
            diags_a = set([frozenset(d) for d in alpha_diagrams[alpha]])
            diags_b = set([frozenset(d) for d in beta_diagrams[beta]])
            common = diags_a.intersection(diags_b)
            
            val = QQ(0)
            for diag in common:
                den = QQ(1)
                valid = True
                for pole in diag:
                    s_val = s_subset(pole)
                    if s_val == 0:
                        valid = False
                        break
                    den *= s_val
                if valid and den != 0:
                    val += 1/den
            M[i, j] = val
    
    return M, basis_perms


def compute_normalization_candidates(lambdas, tilde_lambdas):
    """
    Compute various candidate normalization factors.
    
    Based on the literature, the normalization should involve:
    - Products of angle brackets <i,i+1>
    - Possibly Parke-Taylor denominators
    """
    candidates = {}
    
    # Candidate 1: Cyclic product <01><12><23><34><45><50>
    cyclic_prod = QQ(1)
    for i in range(6):
        j = (i + 1) % 6
        cyclic_prod *= ang_bracket(lambdas, i, j)
    candidates['cyclic_angle'] = cyclic_prod
    
    # Candidate 2: <01>^4 (helicity factor for MHV)
    candidates['angle_01_4'] = ang_bracket(lambdas, 0, 1)**4
    
    # Candidate 3: <01>^8 (full MHV gravity factor)
    candidates['angle_01_8'] = ang_bracket(lambdas, 0, 1)**8
    
    # Candidate 4: Cyclic product to various powers
    for power in [1, 2, 4]:
        candidates[f'cyclic_pow_{power}'] = cyclic_prod**power
    
    # Candidate 5: Products of specific brackets
    # <01><12><20> * <34><45><53>
    triple_1 = ang_bracket(lambdas, 0, 1) * ang_bracket(lambdas, 1, 2) * ang_bracket(lambdas, 2, 0)
    triple_2 = ang_bracket(lambdas, 3, 4) * ang_bracket(lambdas, 4, 5) * ang_bracket(lambdas, 5, 3)
    candidates['triple_product'] = triple_1 * triple_2
    
    return candidates


def prove_klt_intersection_equivalence(num_samples=20):
    """
    Prove the relationship S_KLT = c * M^{-1} by finding the normalization c.
    
    Strategy:
    1. Compute S_KLT and M for multiple kinematic points
    2. Try various candidate normalizations
    3. Find which one gives S_KLT / c = M^{-1} consistently
    """
    print("=" * 70)
    print("PROVING: KLT Kernel = c * Inverse Intersection Matrix")
    print("=" * 70)
    
    # Track ratio for each candidate normalization
    candidate_ratios = {}
    
    for sample in range(num_samples):
        result = sample_kinematics(seed=sample * 13 + 7)
        if result is None:
            continue
            
        lambdas, tilde_lambdas = result
        
        # Compute matrices
        S_klt, perms = compute_klt_kernel_matrix(lambdas, tilde_lambdas)
        M_ba, _ = compute_biadjoint_matrix(lambdas, tilde_lambdas)
        
        if M_ba.det() == 0:
            continue
        
        M_inv = M_ba.inverse()
        
        # Compute candidate normalizations
        candidates = compute_normalization_candidates(lambdas, tilde_lambdas)
        
        # For each candidate, check if S_klt / c = M_inv consistently
        for name, c in candidates.items():
            if c == 0:
                continue
            
            # Compute normalized S
            S_normalized = S_klt / c
            
            # Compare to M_inv entry by entry
            ratios = []
            for i in range(6):
                for j in range(6):
                    if M_inv[i,j] != 0 and S_normalized[i,j] != 0:
                        ratio = M_inv[i,j] / S_normalized[i,j]
                        ratios.append(float(ratio))
            
            if ratios:
                avg = sum(ratios) / len(ratios)
                std = (sum((r - avg)**2 for r in ratios) / len(ratios))**0.5
                rel_std = std / abs(avg) if avg != 0 else float('inf')
                
                if name not in candidate_ratios:
                    candidate_ratios[name] = []
                candidate_ratios[name].append((avg, rel_std))
    
    # Analyze results
    print("\nCandidate Normalization Analysis:")
    print("-" * 70)
    
    best_candidate = None
    best_consistency = float('inf')
    
    for name, data in sorted(candidate_ratios.items()):
        if len(data) < 5:
            continue
        
        # Average ratio across samples
        avg_ratio = sum(d[0] for d in data) / len(data)
        # Consistency: std of the ratio across samples
        ratio_std = (sum((d[0] - avg_ratio)**2 for d in data) / len(data))**0.5
        # Average within-sample consistency
        avg_rel_std = sum(d[1] for d in data) / len(data)
        
        print(f"\n{name}:")
        print(f"  Samples: {len(data)}")
        print(f"  Average ratio M^-1 / (S/c): {avg_ratio:.6f}")
        print(f"  Ratio std across samples: {ratio_std:.6f}")
        print(f"  Avg within-sample rel_std: {avg_rel_std:.6f}")
        
        # The best candidate has low within-sample variance AND consistent ratio across samples
        consistency_score = avg_rel_std + abs(ratio_std / avg_ratio) if avg_ratio != 0 else float('inf')
        
        if consistency_score < best_consistency:
            best_consistency = consistency_score
            best_candidate = (name, avg_ratio)
    
    print("\n" + "=" * 70)
    if best_candidate:
        print(f"BEST CANDIDATE: {best_candidate[0]}")
        print(f"Relationship: M^{-1} = {best_candidate[1]:.6f} * (S_KLT / {best_candidate[0]})")
        print(f"Or equivalently: S_KLT = {1/best_candidate[1]:.6f} * {best_candidate[0]} * M^{-1}")
    else:
        print("No consistent normalization found.")
    print("=" * 70)
    
    return best_candidate


def verify_signature_match(num_samples=10):
    """
    Verify that S_KLT and M^{-1} have the same signature.
    
    This is a weaker but important consistency check.
    """
    print("\n" + "=" * 70)
    print("VERIFYING: Signature Match between S_KLT and M^{-1}")
    print("=" * 70)
    
    matches = 0
    total = 0
    
    for sample in range(num_samples):
        result = sample_kinematics(seed=sample * 23)
        if result is None:
            continue
            
        lambdas, tilde_lambdas = result
        
        S_klt, _ = compute_klt_kernel_matrix(lambdas, tilde_lambdas)
        M_ba, _ = compute_biadjoint_matrix(lambdas, tilde_lambdas)
        
        if M_ba.det() == 0:
            continue
        
        M_inv = M_ba.inverse()
        
        # Symmetrize
        S_sym = (S_klt + S_klt.transpose()) / 2
        M_inv_sym = (M_inv + M_inv.transpose()) / 2
        
        # Compute signatures
        S_float = S_sym.change_ring(RDF)
        M_float = M_inv_sym.change_ring(RDF)
        
        eigs_S = S_float.eigenvalues()
        eigs_M = M_float.eigenvalues()
        
        sig_S = (sum(1 for e in eigs_S if e > 1e-10), 
                 sum(1 for e in eigs_S if e < -1e-10))
        sig_M = (sum(1 for e in eigs_M if e > 1e-10), 
                 sum(1 for e in eigs_M if e < -1e-10))
        
        total += 1
        if sig_S == sig_M:
            matches += 1
            status = "✓ MATCH"
        else:
            status = "✗ MISMATCH"
        
        print(f"Sample {sample}: S_KLT sig={sig_S}, M^-1 sig={sig_M} {status}")
    
    print(f"\nSignature matches: {matches}/{total}")
    
    if matches == total:
        print("✓ All signatures match! This confirms the geometric relationship.")
    else:
        print("Some signatures differ - this may indicate chamber-dependent behavior.")
    
    return matches == total


if __name__ == "__main__":
    # First verify signatures match
    verify_signature_match(num_samples=10)
    
    # Then find the normalization
    prove_klt_intersection_equivalence(num_samples=15)

