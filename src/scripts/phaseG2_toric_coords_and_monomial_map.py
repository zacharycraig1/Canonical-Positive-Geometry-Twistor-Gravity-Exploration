import sys
import os
import json
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.toric import compute_lattice_basis, get_toric_exponents

def phaseG2_toric_mapping():
    print("Phase G2: Computing Toric Coordinates...")
    
    cache_dir = ".dcp_cache/phaseG"
    files = sorted([f for f in os.listdir(cache_dir) if f.endswith(".json") and "roots_012" in f])
    
    for fname in files:
        path = os.path.join(cache_dir, fname)
        with open(path, "r") as f:
            data = json.load(f)
            
        print(f"\nProcessing {fname}...")
        vertices = data["vertices"]
        n = data["n"]
        
        # 1. Compute Basis
        d, B_ambient, v0 = compute_lattice_basis(vertices)
        print(f"  Affine dimension: {d} (Expected {data['dim']})")
        
        # 2. Map to exponents
        exponents = get_toric_exponents(vertices, B_ambient, v0)
        
        # 3. Save Map Info
        # We need to save the basis so we can interpret physics variables later
        map_data = {
            "d": int(d),
            "origin": [int(x) for x in v0],
            "basis_rows": [[int(c) for c in row] for row in B_ambient.rows()],
            "exponents": [[int(c) for c in u] for u in exponents]
        }
        
        out_name = fname.replace("P_", "Toric_")
        out_path = os.path.join(cache_dir, out_name)
        with open(out_path, "w") as f:
            json.dump(map_data, f, indent=2)
            
        print(f"  Saved toric map to {out_name}")
        
        # Sanity check: verify all vertices are unique in exponent space
        unique_u = set(tuple(u) for u in exponents)
        if len(unique_u) != len(vertices):
            print("  WARNING: Exponent collision! Basis might be degenerate?")
        else:
            print("  Mapping bijective on vertices.")

if __name__ == "__main__":
    phaseG2_toric_mapping()



