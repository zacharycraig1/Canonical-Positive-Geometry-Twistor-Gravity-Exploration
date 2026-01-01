import sys
import os
import json
import random
from sage.all import *

sys.path.append(os.getcwd())

def phaseG3_check_facets():
    print("Phase G3: Checking Facet Structures...")
    
    cache_dir = ".dcp_cache/phaseG"
    files = sorted([f for f in os.listdir(cache_dir) if f.startswith("P_") and f.endswith(".json")])
    
    for fname in files:
        path = os.path.join(cache_dir, fname)
        with open(path, "r") as f:
            data = json.load(f)
            
        print(f"\nPolytope {fname}:")
        facets = data["facets_inequalities"]
        dim = data["dim"]
        
        # Analyze facets
        # Inequality: b + a.x >= 0
        # For Newton polytope in R^d_>=0, usually facets are x_i >= 0 plus others.
        
        count_coord = 0
        count_nontrivial = 0
        
        for f in facets:
            b = f[0]
            a = f[1:]
            # Check if it's a coordinate plane (one 1, rest 0)
            non_zeros = [x for x in a if x != 0]
            if len(non_zeros) == 1 and (non_zeros[0] == 1 or non_zeros[0] == -1):
                count_coord += 1
            else:
                count_nontrivial += 1
                
        print(f"  Facets: {len(facets)}")
        print(f"  Coordinate-like: {count_coord}")
        print(f"  Non-trivial: {count_nontrivial}")
        
if __name__ == "__main__":
    phaseG3_check_facets()






