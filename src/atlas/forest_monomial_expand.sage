import sys
import os
import itertools
import json
from sage.all import *

# Path setup
sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents, enumerate_rooted_forests

def expand_forest_monomials():
    print("Expanding Forest Polynomials Symbolically...")
    
    n = 6
    roots_list = list(itertools.combinations(range(n), 3))
    
    # We want to check the "Matrix Tree Theorem" relation.
    # For a set of roots R, the reduced Laplacian determinant det(L_R)
    # should expand to the sum of forests rooted at R.
    
    # We will verify this for all 20 charts.
    # And we will output the list of monomials for each chart.
    # This serves as a "fingerprint" for the chart.
    
    results = []
    
    for idx, roots in enumerate(roots_list):
        print(f"Checking Chart {idx+1}/{len(roots_list)}: Roots {roots}")
        
        # 1. Enumerate Forests explicitly
        # This gives us the "Target" monomials
        forests = enumerate_rooted_forests(n, roots)
        # Each forest is a list of edges (u,v).
        # We can represent monomial as a tuple of sorted edges.
        
        target_monomials = set()
        for f in forests:
            # Sort edges
            edges = tuple(sorted([tuple(sorted(e)) for e in f]))
            target_monomials.add(edges)
            
        print(f"  Explicit Forests: {len(target_monomials)}")
        
        # 2. Symbolic Determinant
        # Variables z_ij
        edge_vars = {}
        var_names = []
        for i in range(n):
            for j in range(i+1, n):
                name = f"z_{i}_{j}"
                var_names.append(name)
                edge_vars[(i,j)] = name
                edge_vars[(j,i)] = name
                
        R = PolynomialRing(QQ, var_names)
        z = R.gens_dict()
        
        # Build Laplacian L (n x n)
        # L_uv = -z_uv
        # L_uu = sum_k z_uk
        
        L = matrix(R, n, n)
        for i in range(n):
            diag = 0
            for j in range(n):
                if i == j: continue
                # edge (i,j)
                if i < j: uv = (i,j)
                else: uv = (j,i)
                val = z[edge_vars[uv]]
                
                L[i,j] = -val
                diag += val
            L[i,i] = diag
            
        # Reduced Laplacian L_R
        # Remove rows/cols corresponding to roots
        non_roots = sorted([i for i in range(n) if i not in roots])
        
        # We need the minor corresponding to removing roots?
        # The Matrix Tree Theorem says:
        # Sum of forests rooted at R = det(L[rows=non_roots, cols=non_roots])?
        # Wait, the theorem usually says:
        # To get trees rooted at r, remove row/col r.
        # To get forests rooted at {r1...rk}, we remove rows/cols {r1...rk}?
        # Let's verify.
        
        indices = non_roots
        L_R = L[indices, indices] # submatrix
        
        det_poly = L_R.det()
        
        # Extract monomials from poly
        det_monomials = set()
        for exponent, coeff in det_poly.dict().items():
            # exponent is tuple of powers for variables in generic order
            # map back to edges
            
            # Reconstruct edge list
            # exponent is tuple corresponding to R.gens()
            # R.gens() are z_0_1, z_0_2 ...
            
            edges = []
            gens = R.gens()
            for k, power in enumerate(exponent):
                if power > 0:
                    # Which edge?
                    # var_names[k] is "z_u_v"
                    name = var_names[k]
                    # parse
                    parts = name.split('_')
                    u, v = int(parts[1]), int(parts[2])
                    for _ in range(power):
                         edges.append((u,v))
            
            # Sort
            edges_tuple = tuple(sorted(edges))
            det_monomials.add(edges_tuple)
            
            # Verify coeff is 1
            if abs(coeff - 1) > 1e-9:
                print(f"  Warning: Monomial {edges_tuple} has coeff {coeff}")
                
        print(f"  Determinant Monomials: {len(det_monomials)}")
        
        # Compare
        if target_monomials == det_monomials:
            print("  MATCH: Determinant matches Forest Sum exactly.")
            status = "MATCH"
        else:
            print(f"  MISMATCH: Target {len(target_monomials)} vs Det {len(det_monomials)}")
            status = "MISMATCH"
            # Difference?
            diff1 = target_monomials - det_monomials
            diff2 = det_monomials - target_monomials
            if diff1: print(f"    Missing in Det: {list(diff1)[:3]}...")
            if diff2: print(f"    Extra in Det: {list(diff2)[:3]}...")
            
        results.append({
            "roots": list(roots),
            "status": status,
            "count": len(det_monomials)
        })
        
    # Save results
    with open("RESULTS/forest_monomial_check.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    expand_forest_monomials()

