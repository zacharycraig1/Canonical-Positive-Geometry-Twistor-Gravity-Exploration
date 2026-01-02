#!/usr/bin/env sage
"""
Explicit Signed Region in Kinematic Space

Goal: Find the explicit set of inequalities that define the "signed region"
whose signed canonical form gives the gravity amplitude.

For positive geometry:
    X = { z : f_1(z) > 0, f_2(z) > 0, ..., f_k(z) > 0 }
    Ω(X) = amplitude

For signed geometry:
    X = { z : various inequalities with assigned signs }
    Ω_signed(X) = amplitude

The key insight is that the signs come from the kinematic weights w_ij.
The signed region should be defined by conditions on these weights.
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


def analyze_weight_sign_patterns(num_samples=50):
    """
    Analyze the patterns of w_ij signs and how they relate to forest signs.
    
    The weights w_ij = [ij]/⟨ij⟩ have signs that vary with kinematics.
    We want to understand how these signs determine the 54/54 split.
    """
    print("=" * 70)
    print("WEIGHT SIGN PATTERN ANALYSIS")
    print("=" * 70)
    
    n = 6
    
    # Track which w_ij are positive vs negative
    sign_patterns = []
    forest_splits = []
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 59)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        # Compute w_ij signs
        w_signs = {}
        for i in range(n):
            for j in range(i+1, n):
                ang = ang_bracket(lambdas, i, j)
                sq = sq_bracket(tilde_lambdas, i, j)
                if ang != 0:
                    w = sq / ang
                    w_signs[(i, j)] = 1 if w > 0 else -1
        
        # Convert to pattern
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        pattern = tuple(w_signs.get(e, 0) for e in edges)
        sign_patterns.append(pattern)
        
        # Count positive w_ij
        n_pos = sum(1 for s in pattern if s > 0)
        n_neg = sum(1 for s in pattern if s < 0)
        forest_splits.append((n_pos, n_neg))
    
    # Analyze patterns
    print(f"\nSamples: {len(sign_patterns)}")
    print(f"Total edges: {len(edges)} = C(6,2) = 15")
    
    # Distribution of positive weights
    print("\nDistribution of positive w_ij:")
    from collections import Counter
    pos_counts = Counter(s[0] for s in forest_splits)
    for n_pos, count in sorted(pos_counts.items()):
        pct = 100.0 * count / len(forest_splits)
        print(f"  {n_pos} positive: {count} ({pct:.1f}%)")
    
    # Are certain edges more likely to be positive?
    print("\nPer-edge positivity rates:")
    edge_pos_counts = {e: 0 for e in edges}
    for pattern in sign_patterns:
        for i, e in enumerate(edges):
            if pattern[i] > 0:
                edge_pos_counts[e] += 1
    
    for e in edges:
        rate = float(edge_pos_counts[e]) / len(sign_patterns)
        print(f"  w_{e[0]}{e[1]}: {rate:.1%} positive")
    
    return sign_patterns


def define_signed_region():
    """
    Define the signed region conceptually.
    
    The signed region for gravity is NOT a simple positive region.
    Instead, it's the entire kinematic space with a signed measure.
    
    Conceptually:
        Ω_signed = Σ_chambers ε(chamber) × Ω(chamber)
        
    Where chambers are defined by signs of w_ij.
    """
    print("\n" + "=" * 70)
    print("SIGNED REGION DEFINITION")
    print("=" * 70)
    
    print("""
For gravity, the "signed region" is the full kinematic space,
but with a SIGNED MEASURE determined by the edge weights w_ij.

CHAMBERS:
    The kinematic space decomposes into 2^15 chambers based on
    the signs of the 15 weights w_ij = [ij]/⟨ij⟩.
    
SIGNED CANONICAL FORM:
    In each chamber, the canonical form is:
        Ω_chamber = ε(chamber) × |Ω|
    
    where ε(chamber) depends on which w_ij are positive.
    
AMPLITUDE:
    M = ∫_X Ω_signed = Σ_chambers Σ_forests-in-chamber ε(F) × ω(F)
    
    But since forests contribute regardless of chamber, the sum
    is just:
        M = Σ_forests ε(F) × ω(F)
    
    where ε(F) is determined by the sign rule.

KEY INSIGHT:
    Unlike positive geometry where we need X to be a polytope,
    for signed geometry we work on the full space but track signs.
    
    The "geometry" is in the SIGN STRUCTURE, not the region.
""")


def analyze_chamber_to_forest_connection():
    """
    Analyze how chambers (w_ij sign patterns) connect to forest signs.
    
    For each chamber, how many forests have + vs - sign?
    """
    print("\n" + "=" * 70)
    print("CHAMBER TO FOREST SIGN CONNECTION")
    print("=" * 70)
    
    load('src/signed_geometry/forest_sign_rule.sage')
    
    n = 6
    roots = (0, 1, 2)
    all_forests = list(enumerate_rooted_forests(n, roots))
    
    # For each forest, determine which w_ij it uses
    forest_edges = {}
    for forest in all_forests:
        forest_edges[forest] = set(forest)
    
    print(f"Total forests: {len(all_forests)}")
    print(f"Edges per forest: 3")
    
    # For a random kinematic sample, analyze
    result = sample_spinor_helicity_conserving(n=6, seed=42)
    if result is None:
        return
    
    lambdas, tilde_lambdas = result
    x_spinor = vector(QQ, [1, 2])
    y_spinor = vector(QQ, [3, 1])
    
    # Get w_ij signs
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
    
    # For each forest, compute its sign from our sign rule
    # and verify it matches the actual sign
    forest_signs = {}
    
    for forest in all_forests:
        edges = list(forest)
        
        # Sign from (-1)^|E|
        sign_from_edges = (-1)**len(edges)
        
        # Sign from w product
        sign_from_w = 1
        for u, v in edges:
            if w[(u, v)] < 0:
                sign_from_w *= -1
        
        # Sign from C product
        degree = {i: 0 for i in range(n)}
        for u, v in edges:
            degree[u] += 1
            degree[v] += 1
        
        sign_from_C = 1
        for v in range(n):
            if degree[v] % 2 == 1 and C[v] < 0:
                sign_from_C *= -1
        
        # Total sign
        forest_signs[forest] = sign_from_edges * sign_from_w * sign_from_C
    
    pos_forests = sum(1 for s in forest_signs.values() if s > 0)
    neg_forests = sum(1 for s in forest_signs.values() if s < 0)
    
    print(f"\nForest sign distribution: {pos_forests}+, {neg_forests}-")
    
    # Group forests by which edges they use
    # and analyze if certain edge combinations lead to certain signs
    print("\nForests grouped by edge pattern:")
    
    # Count how many times each edge appears in positive vs negative forests
    edge_in_pos = {e: 0 for e in w.keys()}
    edge_in_neg = {e: 0 for e in w.keys()}
    
    for forest, sign in forest_signs.items():
        for e in forest:
            if sign > 0:
                edge_in_pos[e] += 1
            else:
                edge_in_neg[e] += 1
    
    print("\nEdge appearance in + vs - forests:")
    for e in sorted(w.keys()):
        w_sign = "+" if w[e] > 0 else "-"
        print(f"  w_{e[0]}{e[1]} ({w_sign}): {edge_in_pos[e]}+ / {edge_in_neg[e]}-")


if __name__ == "__main__":
    analyze_weight_sign_patterns(num_samples=30)
    define_signed_region()
    analyze_chamber_to_forest_connection()

