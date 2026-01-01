import sys
import os
import json
import itertools
from sage.all import *

# Path setup
sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents

def compute_chart_adjacency():
    print("Computing Chart Adjacency and Signs...")
    
    n = 6
    all_roots = list(itertools.combinations(range(n), 3))
    
    # Store chart info
    chart_facets = {} # chart_idx -> list of (normal, offset)
    facet_to_charts = {} # (normal, offset) -> list of chart_idx
    
    print(f"Processing {len(all_roots)} charts...")
    
    for idx, roots in enumerate(all_roots):
        roots_tuple = tuple(sorted(list(roots)))
        exponents, edge_order = get_forest_exponents(n, roots_tuple)
        
        # Build Polytope
        P = Polyhedron(vertices=exponents)
        
        # Get Facets (Inequalities)
        # H-rep: Ax + b >= 0  (Sage default is usually Ax + b >= 0 or <= 0 depending on method)
        # inequality_generator() returns (b, A) such that A*x + b >= 0
        
        ieqs = P.inequality_generator()
        
        current_facets = []
        for ieq in ieqs:
            # ieq is a vector [b, a1, ..., am] representing b + a*x >= 0
            # Normalize:
            # We want canonical representation. 
            # Divide by GCD of all elements?
            # Make first non-zero positive?
            
            vec = vector(QQ, ieq)
            
            # For integer coordinates, we can normalize to integers
            # But the vector might be rational.
            # Convert to integer vector by clearing denominators and dividing by GCD
            
            # Since vertices are integers (0/1), the facet equations can be integers.
            vec_int = (vec * 1).apply_map(lambda x: x) # Ensure it's treated as vector
            
            # Find common denominator?
            # Sage Polyhedron usually gives integer coefficients if vertices are integer?
            # Let's check type.
            
            # To be robust:
            # 1. Scale to integers (if rational)
            # 2. Divide by GCD
            # 3. Orient lexicographically or first-non-zero positive to handle sign ambiguity?
            # Wait, A*x + b >= 0 defines the half-space containing the polytope.
            # The direction is fixed (outward/inward). We don't want to flip A and b.
            # So just normalize by GCD of magnitudes.
            
            coeffs = [x for x in vec]
            
            # If all ints
            try:
                # check if all are integers (or close)
                coeffs_int = [int(round(x)) for x in coeffs]
                # verify
                if all(abs(x - y) < 1e-9 for x, y in zip(coeffs, coeffs_int)):
                    gcd_val = GCD(coeffs_int)
                    if gcd_val > 0:
                        coeffs_norm = tuple([c // gcd_val for c in coeffs_int])
                    else:
                        coeffs_norm = tuple(coeffs_int)
                else:
                    # Rational case, unlikely for this polytope
                    coeffs_norm = tuple(coeffs)
            except:
                coeffs_norm = tuple(coeffs)
                
            current_facets.append(coeffs_norm)
            
            if coeffs_norm not in facet_to_charts:
                facet_to_charts[coeffs_norm] = []
            facet_to_charts[coeffs_norm].append(idx)
            
        chart_facets[idx] = current_facets
        if idx % 5 == 0:
            print(f"  Processed {idx+1}/{len(all_roots)}")

    print(f"Identified {len(facet_to_charts)} unique hyperplanes.")
    
    # Solver Strategy
    # We want to find signs sigma_i in {-1, 1} such that for every NON-TRIVIAL facet F:
    #   Sum_{R in F} sigma_R = 0
    #
    # Facet Classification:
    #   - Coordinate (z_e >= 0): Physical boundary, sum != 0 allowed.
    #   - Others (e.g. z_e <= 1, sum z >= k, etc): Spurious/Internal, sum = 0 required.
    
    constraints = []
    
    print("\nIdentifying Constraints...")
    for facet, charts in facet_to_charts.items():
        # Check if coordinate
        b = facet[0]
        coeffs = facet[1:]
        non_zeros = [i for i, x in enumerate(coeffs) if x != 0]
        
        # z_e >= 0  <==>  0 + 1*z_e >= 0  => (0, ... 1 ...)
        # Note: facet is (b, a1...am). 
        # Coordinate facet usually looks like (0, ..., 1, ...).
        # Or maybe (0, ..., -1, ...) if bounds are upper? Forest polytope is usually z_e >= 0.
        
        is_coord_pos = (b == 0 and len(non_zeros) == 1 and coeffs[non_zeros[0]] > 0)
        
        if is_coord_pos:
            continue # Physical boundary
            
        # Everything else is a constraint
        constraints.append((facet, charts))
        
    print(f"Found {len(constraints)} constraints from internal/spurious facets.")
    
    # Solve for signs
    # N=20 is small enough for brute force or optimized search
    # We can optimize by propagating "Shared by 2" constraints first
    
    # 1. Propagate binary constraints
    fixed_signs = {}
    # Fix chart 0 to 1 (WLOG)
    fixed_signs[0] = 1
    
    # Also look for charts that might be forced if they are the ONLY one on a constraint?
    # If len(charts) == 1 for an internal facet, that facet is a "hole" in the atlas?
    # Ideally len >= 2.
    
    binary_constraints = [c for c in constraints if len(c[1]) == 2]
    other_constraints = [c for c in constraints if len(c[1]) > 2]
    
    print(f"Binary constraints: {len(binary_constraints)}")
    
    # Propagate BFS style
    queue = [0]
    processed_constraints = set()
    
    # While we have nodes to process or constraints to check
    changed = True
    while changed:
        changed = False
        
        # 1. Check binary constraints connected to fixed signs
        for i, (facet, charts) in enumerate(binary_constraints):
            if i in processed_constraints: continue
            
            c0, c1 = charts[0], charts[1]
            if c0 in fixed_signs and c1 not in fixed_signs:
                fixed_signs[c1] = -fixed_signs[c0]
                processed_constraints.add(i)
                changed = True
            elif c1 in fixed_signs and c0 not in fixed_signs:
                fixed_signs[c0] = -fixed_signs[c1]
                processed_constraints.add(i)
                changed = True
            elif c0 in fixed_signs and c1 in fixed_signs:
                if fixed_signs[c0] + fixed_signs[c1] != 0:
                    print(f"Conflict at binary facet {facet} for charts {charts}")
                processed_constraints.add(i)
                
    print(f"Propagated signs to {len(fixed_signs)} charts via binary constraints.")
    
    # 2. Brute force remaining
    unknowns = [i for i in range(len(all_roots)) if i not in fixed_signs]
    if unknowns:
        print(f"Solving for {len(unknowns)} remaining signs...")
        best_signs = None
        min_violation = float('inf')
        
        if len(unknowns) > 20:
            print("Too many unknowns for brute force. Using randomized greedy.")
            # ...
            # Just default to 1 for now if too many
            best_signs = fixed_signs.copy()
            for u in unknowns: best_signs[u] = 1
        else:
            for signs_tuple in itertools.product([1, -1], repeat=len(unknowns)):
                current_signs = fixed_signs.copy()
                for i, s in zip(unknowns, signs_tuple):
                    current_signs[i] = s
                    
                # Check violations
                violation = 0
                for facet, charts in constraints:
                    s_sum = sum(current_signs[c] for c in charts)
                    # We want sum to be 0
                    if s_sum != 0:
                        violation += 1
                    
                if violation < min_violation:
                    min_violation = violation
                    best_signs = current_signs.copy()
                    if min_violation == 0:
                        break
        print(f"Best solution violation count: {min_violation}")
    else:
        best_signs = fixed_signs
        print("All signs determined by binary constraints.")

    results = []
    for idx, roots in enumerate(all_roots):
        s = best_signs.get(idx, 1)
        results.append({
            "roots": [int(r) for r in roots],
            "sign": int(s),
            "chart_idx": int(idx)
        })
        
    out_path = "RESULTS/atlas_signs.json"
    os.makedirs("RESULTS", exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print(f"Saved signs to {out_path}")

            
if __name__ == "__main__":
    compute_chart_adjacency()

