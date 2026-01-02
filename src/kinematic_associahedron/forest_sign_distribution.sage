# forest_sign_distribution.sage
"""
DETAILED SIGN DISTRIBUTION ANALYSIS

Key question: For a typical kinematic point, how are the 108 forest signs distributed?
- If it's 54/54 (half positive, half negative), there might be a fundamental obstruction
- If it's 107/1 (almost all same sign), we're close and might need a small adjustment

Also: what's the structure of the sign pattern?
"""

from sage.all import *
from itertools import combinations, product
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics

print("="*70)
print("FOREST SIGN DISTRIBUTION ANALYSIS")
print("="*70)

def enumerate_forests(n, roots):
    """Enumerate all spanning forests of K_n with trees rooted at 'roots'."""
    non_roots = [v for v in range(n) if v not in roots]
    forests = []
    
    all_vertices = list(range(n))
    
    for parents in product(all_vertices, repeat=len(non_roots)):
        valid = True
        edges = []
        
        for i, nr in enumerate(non_roots):
            p = parents[i]
            if p == nr:
                valid = False
                break
            edges.append((min(nr, p), max(nr, p)))
        
        if not valid:
            continue
        
        parent_uf = {v: v for v in range(n)}
        
        def find(x):
            if parent_uf[x] != x:
                parent_uf[x] = find(parent_uf[x])
            return parent_uf[x]
        
        def union(x, y):
            px, py = find(x), find(y)
            if px == py:
                return False
            parent_uf[px] = py
            return True
        
        is_forest = True
        for (a, b) in edges:
            if not union(a, b):
                is_forest = False
                break
        
        if not is_forest:
            continue
        
        root_components = [find(r) for r in roots]
        if len(set(root_components)) != len(roots):
            continue
        
        forests.append(tuple(sorted(edges)))
    
    return list(set(forests))


forests = enumerate_forests(6, (0, 1, 2))
print(f"Enumerated {len(forests)} forests\n")

# Collect all edges
all_edges = set()
for F in forests:
    for e in F:
        all_edges.add(e)
all_edges = sorted(all_edges)
print(f"Edges: {all_edges}")

# Test distribution for several kinematic points
print("\n" + "="*70)
print("SIGN DISTRIBUTION FOR SAMPLE KINEMATICS")
print("="*70)

distributions = []

for seed in [42, 123, 456, 789, 2024]:
    try:
        kin = SpinorKinematics.random_rational(n=6, seed=seed)
        lambdas = kin.lambdas
        lambdas_tilde = kin.tilde_lambdas
        
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
        
        if any(v is None for v in w.values()):
            continue
        
        # Compute sign of each forest term
        positive_count = 0
        negative_count = 0
        
        for F in forests:
            prod = 1
            for e in F:
                prod *= w[e]
            if prod > 0:
                positive_count += 1
            else:
                negative_count += 1
        
        ratio = positive_count / len(forests) if len(forests) > 0 else 0
        distributions.append(positive_count)
        
        print(f"\nSeed {seed}:")
        print(f"  Positive: {positive_count}, Negative: {negative_count}")
        print(f"  Ratio: {100*ratio:.1f}%")
        
    except Exception as e:
        print(f"Seed {seed}: Error - {e}")
        continue

print("\n" + "="*70)
print("ANALYSIS")
print("="*70)

if distributions:
    avg = float(sum(distributions)) / len(distributions)
    print(f"\nAverage positive forests: {avg:.1f} out of 108")
    print(f"Range: {min(distributions)} to {max(distributions)}")
    
    # The key insight
    print(f"""
INTERPRETATION:
If the distribution is roughly 50/50 (54 positive, 54 negative), there's 
no simple positive region where the amplitude has definite sign.

If the distribution is skewed (e.g., 80/28), then there might be a region 
where ALL forests have the same sign, just hard to find by random sampling.

From our data: average = {avg:.1f} forests positive out of 108.
""")

# Check: what sign patterns occur for edges?
print("\n" + "="*70)
print("EDGE SIGN PATTERNS")
print("="*70)

# For each sample, what's the sign pattern of the 12 edges?
print("\nAnalyzing which edge sign patterns lead to more consistent forest signs...")

best_consistency = 0
best_seed = None
best_pattern = None

for seed in range(1, 501):
    try:
        kin = SpinorKinematics.random_rational(n=6, seed=seed)
        lambdas = kin.lambdas
        lambdas_tilde = kin.tilde_lambdas
        
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
        
        if any(v is None for v in w.values()):
            continue
        
        # Edge sign pattern
        pattern = tuple('+' if w[e] > 0 else '-' for e in all_edges)
        
        # Count consistent forests
        positive_count = 0
        negative_count = 0
        
        for F in forests:
            prod = 1
            for e in F:
                prod *= w[e]
            if prod > 0:
                positive_count += 1
            else:
                negative_count += 1
        
        consistency = max(positive_count, negative_count)
        
        if consistency > best_consistency:
            best_consistency = consistency
            best_seed = seed
            best_pattern = pattern
            
    except:
        continue

print(f"\nBest consistency found: {best_consistency}/108 forests with same sign")
print(f"Seed: {best_seed}")
print(f"Edge sign pattern: {best_pattern}")

if best_consistency == 108:
    print("\n✅ FOUND CONSISTENT KINEMATICS!")
elif best_consistency >= 100:
    print(f"\n🔶 Very close! {best_consistency}/108 = {best_consistency/108:.1%} consistent")
else:
    print(f"\n⚠️ Best achievable: {best_consistency}/108 = {best_consistency/108:.1%}")
    print("   A fully consistent positive region may not exist for generic kinematics.")

