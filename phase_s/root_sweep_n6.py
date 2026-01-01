
import sys
import os
import json
import itertools
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from phase_s.analysis_utils import analyze_limit_with_retry

def run_root_sweep():
    n = 6
    epsilon = QQ(1)/1000
    
    # All 20 root choices
    all_roots = list(itertools.combinations(range(n), 3))
    
    results = {}
    
    limits_to_check = [
        ("Soft(5)", "soft", {"s_idx": 5}),
        ("Col(3||4)", "collinear", {"i": 3, "j": 4}),
        ("Col(4||5)", "collinear", {"i": 4, "j": 5})
    ]
    
    print(f"Starting sweep over {len(all_roots)} root sets...")
    
    for roots in all_roots:
        roots_tuple = tuple(sorted(list(roots)))
        roots_str = str(roots_tuple)
        print(f"Roots: {roots_str}")
        
        results[roots_str] = {}
        
        for name, ltype, kwargs in limits_to_check:
            res = analyze_limit_with_retry(n, list(roots), name, ltype, epsilon, **kwargs)
            
            # Store top active facet (slack < 0.01)
            active = []
            if res:
                for f in res:
                    if f["slack"] < 0.01:
                        active.append(f["description"])
            
            results[roots_str][name] = active
            
    # Save results
    with open("phase_s/root_sweep_results.json", "w") as f:
        json.dump(results, f, indent=2)
        
    # Analyze coverage
    # We want to see if for each limit, there EXISTS a root set where it hits a facet.
    # And specifically for Col(4||5), does it hit z_4_5 >= 0?
    
    print("\nCoverage Analysis:")
    
    for name, _, _ in limits_to_check:
        print(f"\nLimit {name}:")
        hit_sets = []
        descriptions = set()
        
        for roots, data in results.items():
            if data[name]:
                hit_sets.append(roots)
                for d in data[name]:
                    descriptions.add(d)
                    
        print(f"  Captured by {len(hit_sets)}/20 root sets")
        print(f"  Active Facets types: {list(descriptions)}")
        
        # Check specifically for Col(4||5)
        if name == "Col(4||5)":
            if any("z_4_5" in d for d in descriptions):
                print("  -> FOUND z_4_5 active!")
            else:
                print("  -> z_4_5 NOT found active.")

if __name__ == "__main__":
    run_root_sweep()



