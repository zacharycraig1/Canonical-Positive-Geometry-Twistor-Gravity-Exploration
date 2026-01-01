import json
import os
import sys
from collections import defaultdict

def analyze_facets_n6():
    print("Analyzing Facets of n=6 Forest Polytope...")
    
    cache_file = ".dcp_cache/phaseG/P_n6_roots_012.json"
    if not os.path.exists(cache_file):
        print("Error: Cache file not found. Run phaseG1_build_forest_polytope_data.py first.")
        return
        
    with open(cache_file, "r") as f:
        data = json.load(f)
        
    edges = data['edges_ordered'] # List of [u, v]
    inequalities = data['facets_inequalities']
    
    print(f"Loaded {len(inequalities)} facets for n={data['n']} roots={data['roots']}.")
    
    # Format: [b, a1, a2, ... aN] where N is number of edges (15 for n=6)
    # Inequality: b + sum(ai * xi) >= 0
    
    # Classify facets
    # Types:
    # 1. Coordinate hyperplanes: xi >= 0  (b=0, one ai=1, others 0)
    # 2. Upper bounds: xi <= 1 (b=1, one ai=-1) -- unlikely for Newton polytopes usually, but possible
    # 3. Sum constraints: sum xi <= k or sum xi >= k
    
    facet_types = defaultdict(list)
    
    for idx, ieq in enumerate(inequalities):
        b = ieq[0]
        coeffs = ieq[1:]
        
        # Check for simple coordinate facets
        non_zero = [(i, c) for i, c in enumerate(coeffs) if c != 0]
        
        if len(non_zero) == 1:
            idx_edge, coeff = non_zero[0]
            u, v = edges[idx_edge]
            edge_str = f"z_{u}_{v}"
            
            if b == 0 and coeff == 1:
                desc = f"{edge_str} >= 0"
                facet_types["Non-negativity"].append(desc)
            elif b == 0 and coeff == -1:
                desc = f"{edge_str} <= 0"
                facet_types["Non-positivity"].append(desc)
            elif coeff == -1:
                desc = f"{edge_str} <= {b}"
                facet_types["Upper Bound"].append(desc)
            elif coeff == 1:
                desc = f"{edge_str} >= {-b}"
                facet_types["Lower Bound"].append(desc)
            else:
                desc = f"{coeff}*{edge_str} + {b} >= 0"
                facet_types["Scaled Bound"].append(desc)
        else:
            # Complex facet
            # Construct string representation
            terms = []
            for i, c in non_zero:
                u, v = edges[i]
                if c == 1:
                    terms.append(f"z_{u}_{v}")
                elif c == -1:
                    terms.append(f"-z_{u}_{v}")
                else:
                    terms.append(f"{c}z_{u}_{v}")
            
            # Form: sum terms >= -b
            lhs = " + ".join(terms).replace("+ -", "- ")
            desc = f"{lhs} >= {-b}"
            facet_types["Complex"].append(desc)
            
    # Report
    print("\n--- Facet Classification ---")
    for category, facets in facet_types.items():
        print(f"\n{category} ({len(facets)}):")
        for f in sorted(facets):
            print(f"  {f}")
            
    # Specific Check: Do we see subtour elimination or similar?
    # Forest polytope facets usually correspond to:
    # x_e >= 0
    # sum_{e in S} x_e <= r(S) for closed sets S?
    # Or sum_{e in S} x_e >= ...?
    
    # Interpretation for Physics
    # z_ij -> 0 corresponds to non-negativity z_ij >= 0 facets.
    # Are there others?
    
    print("\n--- Summary ---")
    print(f"Total Facets: {len(inequalities)}")
    print(f"Coordinate Facets: {len(facet_types['Non-negativity'])}")
    print(f"Other Facets: {len(inequalities) - len(facet_types['Non-negativity'])}")

if __name__ == "__main__":
    analyze_facets_n6()




