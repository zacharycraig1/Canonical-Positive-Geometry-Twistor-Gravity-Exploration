# analyze_atlas_structure.sage
"""
DEEP ANALYSIS: Atlas Structure for Gravity Positive Geometry

The key observation from the atlas scan:
- Each chart with roots R = {r1, r2, r3} has poles at s_ij where i,j ∉ R
- For n=6, each chart has poles at C(3,2) = 3 pairs
- Total charts = C(6,3) = 20
- Total poles = C(6,2) = 15

Critical question: What are the atlas coefficients c_R such that
    M_gravity = Σ_R c_R × Ω_R ?

Strategy:
1. Build the "incidence matrix" A where A[pole][chart] = 1 if chart captures pole
2. The residue at each pole must equal the gravity residue
3. This gives 15 constraints on 20 unknowns -> underdetermined
4. Look for patterns (symmetry, BCJ-like structure)
"""

from sage.all import *
from itertools import combinations
import json

# All possible roots (charts)
n = 6
all_charts = list(combinations(range(n), 3))
print(f"Total charts: {len(all_charts)}")

# All poles s_ij
all_poles = [(i, j) for i in range(n) for j in range(i+1, n)]
print(f"Total poles: {len(all_poles)}")

# Build incidence matrix: which charts capture which poles?
# A chart R captures pole s_ij if i ∉ R AND j ∉ R
incidence = {}
for R in all_charts:
    R_set = set(R)
    captured = []
    for (i, j) in all_poles:
        if i not in R_set and j not in R_set:
            captured.append((i, j))
    incidence[R] = captured
    print(f"Chart {R} captures poles: {captured}")

print("\n" + "="*60)
print("INCIDENCE MATRIX")
print("="*60)

# Build matrix A: rows = poles, cols = charts
A = matrix(QQ, len(all_poles), len(all_charts))
for col_idx, R in enumerate(all_charts):
    R_set = set(R)
    for row_idx, (i, j) in enumerate(all_poles):
        if i not in R_set and j not in R_set:
            A[row_idx, col_idx] = 1

print("\nMatrix dimensions:", A.dimensions())
print("Rank of incidence matrix:", A.rank())

# Each pole is captured by how many charts?
print("\n" + "="*60)
print("POLE COVERAGE")
print("="*60)
for row_idx, pole in enumerate(all_poles):
    count = sum(A[row_idx, col_idx] for col_idx in range(len(all_charts)))
    print(f"Pole s_{pole[0]}{pole[1]} covered by {count} charts")

# Find the null space of A^T
# Solutions are vectors v such that A^T v = 0
# This tells us constraints on the chart coefficients
print("\n" + "="*60)
print("ANALYZING SOLUTION SPACE")
print("="*60)

# For the atlas to work: each pole must have residue from charts = residue of gravity
# If gravity has residue r_p at pole p, and chart R has residue 1 at pole p (when R captures p),
# then we need: Σ_R c_R × (1 if R captures p else 0) = r_p
#
# But the issue is: gravity residue at s_ij is the product of lower-point amplitudes,
# which depends on the kinematics, not just the pole structure.
#
# The simplest hypothesis: all residues are the same (with appropriate normalization),
# so c_R = 1/k where k = number of charts capturing each pole

# Check if all poles have same coverage
coverage = [sum(A[row_idx, col_idx] for col_idx in range(len(all_charts))) for row_idx in range(len(all_poles))]
print(f"Coverage distribution: {set(coverage)}")
if len(set(coverage)) == 1:
    k = coverage[0]
    print(f"All poles covered by exactly {k} charts")
    print(f"Naive atlas coefficients: c_R = 1/{k} = {1/k}")
else:
    print("Non-uniform coverage - need more sophisticated analysis")

# Check S6 symmetry
print("\n" + "="*60)
print("SYMMETRY ANALYSIS")
print("="*60)

# The 6-point amplitude should have S6 symmetry (relabeling particles)
# For the atlas to respect this, all charts should have the same coefficient
# This is consistent with c_R = 1/k for all R

# But wait - the actual test shows omega_slope = -1 for matching charts
# and gravity_slope = -1. This suggests c_R = 1, but then we overcount!

# Let's check: if c_R = 1 for all charts, what's the total contribution?
# At each pole, we get k copies of the chart form
# So the total residue would be k × (correct residue)
# To get the correct amplitude, we need c_R = 1/k

# For n=6 with 3 roots, the complement has 3 elements, giving C(3,2) = 3 poles per chart
# And each pole is in C(n-2, 3-2) = C(4, 1) = 4 charts' complements
# Wait, let me recalculate...

# Chart R = {r1, r2, r3}
# Complement = {c1, c2, c3}
# Poles captured = {(c1,c2), (c1,c3), (c2,c3)} = 3 poles

# For pole (i,j), which charts capture it?
# Charts R where i ∉ R and j ∉ R
# = all subsets of {0,...,5} \ {i,j} of size 3
# = C(4, 3) = 4 charts

print(f"\nEach pole is captured by C(4,3) = 4 charts")
print(f"Therefore, naive atlas coefficient should be c_R = 1/4 for all R")

# Now let's verify: do all 20 charts with c_R = 1/4 give the right amplitude?
# This requires checking against actual gravity amplitude values

print("\n" + "="*60)
print("MINIMAL COVER ANALYSIS")
print("="*60)

# The 7-chart minimal cover: each pole captured exactly once
minimal_cover = [
    (0, 1, 2),  # covers s_34, s_35, s_45
    (0, 3, 5),  # covers s_12, s_24, s_14
    (1, 3, 4),  # covers s_02, s_25, s_05
    (2, 4, 5),  # covers s_03, s_01, s_13
    (0, 1, 4),  # covers s_23
    (0, 2, 3),  # covers s_15
    (1, 2, 5),  # covers s_04
]

# Verify coverage
covered_poles = set()
for R in minimal_cover:
    R_set = set(R)
    for (i, j) in all_poles:
        if i not in R_set and j not in R_set:
            if (i, j) in covered_poles:
                print(f"WARNING: Pole ({i},{j}) covered multiple times!")
            covered_poles.add((i, j))

if len(covered_poles) == len(all_poles):
    print(f"Minimal cover with {len(minimal_cover)} charts covers all {len(all_poles)} poles exactly once")
else:
    print(f"Coverage incomplete: {len(covered_poles)}/{len(all_poles)} poles")

# For minimal cover, the coefficients could all be 1 (no overcounting)
print(f"\nMinimal cover hypothesis: M_gravity = Σ_{len(minimal_cover)} Ω_R (coefficient = 1 each)")

print("\n" + "="*60)
print("KEY INSIGHT: TWO POSSIBLE ATLAS FORMS")
print("="*60)

print("""
OPTION 1: Full symmetric sum (20 charts)
    M_gravity = (1/4) × Σ_{all 20 charts} Ω_R
    
    This respects S6 symmetry manifestly.
    Each pole is captured by 4 charts, with coefficient 1/4 each -> residue = 1.
    
OPTION 2: Minimal cover (7 charts)
    M_gravity = Σ_{7 minimal cover charts} Ω_R
    
    Each pole captured exactly once with coefficient 1.
    But this breaks manifest S6 symmetry!
    
The physics should prefer OPTION 1 for symmetry reasons.
But the actual test needs numerical verification.
""")

# Write the pole-chart mapping for later use
pole_chart_map = {}
for row_idx, pole in enumerate(all_poles):
    pole_str = f"s_{pole[0]}{pole[1]}"
    charts_covering = []
    for col_idx, R in enumerate(all_charts):
        if A[row_idx, col_idx] == 1:
            charts_covering.append(R)
    pole_chart_map[pole_str] = charts_covering

with open('results/pole_chart_mapping.json', 'w') as f:
    # Convert tuples to lists for JSON
    json_map = {k: [list(c) for c in v] for k, v in pole_chart_map.items()}
    json.dump(json_map, f, indent=2)

print("\nWrote pole-chart mapping to results/pole_chart_mapping.json")

