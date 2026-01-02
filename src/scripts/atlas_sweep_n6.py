import sys
import os
import itertools
import json
from sage.all import Polyhedron, QQ

sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.facet_to_subset import generate_subset_vectors, normalize_vector

def atlas_sweep():
    n = 6
    all_roots = list(itertools.combinations(range(n), 3))
    print(f"Sweeping atlas for n={n} over {len(all_roots)} root sets...")
    
    seen_subsets = set()
    
    for idx, roots in enumerate(all_roots):
        roots = list(roots)
        print(f"[{idx+1}/{len(all_roots)}] Processing roots {roots}...")
        
        # 1. Build Polytope
        exponents, edges_ordered = get_forest_exponents(n, roots)
        P = Polyhedron(vertices=exponents, base_ring=QQ)
        
        # 2. Get Facets (Inequalities)
        # Note: P.inequalities() might return many if ambient dim is high
        ineqs = P.inequalities()
        
        # 3. Match to Subsets
        # Generate candidates for THIS edge ordering (depends on roots? No, edges are canonical usually)
        # Check forest_polytope.py: edges_ordered is always canonical (0,1)...(n-2,n-1)
        # Yes.
        candidates = generate_subset_vectors(n, edges_ordered)
        
        for h in ineqs:
            vec = list(h.vector())
            # [b, A...]
            A = vec[1:]
            
            # Normalize A
            norm_A, _ = normalize_vector([int(x) for x in A])
            neg_norm_A, _ = normalize_vector([-int(x) for x in A])
            
            match = None
            if norm_A in candidates:
                match = candidates[norm_A]
            elif neg_norm_A in candidates:
                match = candidates[neg_norm_A]
                
            if match:
                s_tuple = tuple(sorted(match["subset"]))
                seen_subsets.add(s_tuple)
                
    print(f"\nTotal unique subsets seen: {len(seen_subsets)}")
    
    # Expected: All non-trivial cuts?
    # Total subsets of size 2, 3, 4?
    # Size 2: 15
    # Size 3: 20 (halved? 10 pairs? No, subsets)
    # Total non-trivial: 2^6 - 2 - 6 = 56?
    # Let's save the list.
    
    sorted_subsets = sorted(list(seen_subsets), key=lambda x: (len(x), x))
    print("Seen subsets:")
    for s in sorted_subsets:
        print(f"  {s}")
        
    results = {
        "n": n,
        "total_charts": len(all_roots),
        "unique_subsets_count": len(seen_subsets),
        "subsets": [list(s) for s in sorted_subsets]
    }
    
    with open("RESULTS/atlas_sweep_n6.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("Saved RESULTS/atlas_sweep_n6.json")

if __name__ == "__main__":
    atlas_sweep()




