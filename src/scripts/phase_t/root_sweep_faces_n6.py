import itertools
import sys
import os
from sage.all import *

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

from posgeom.forest_polytope import get_forest_exponents

def analyze_root_sweep():
    n = 6
    k = 3 # 3-forests
    all_roots = list(itertools.combinations(range(n), k))
    
    report_lines = []
    report_lines.append("# Root Patch Atlas for n=6")
    report_lines.append(f"Analyzing all {len(all_roots)} root sets for n={n}, k={k}.")
    report_lines.append("")
    
    # Track coverage of edges
    edge_coverage = {}
    
    for roots in all_roots:
        roots_list = list(roots)
        print(f"Analyzing roots: {roots_list}")
        
        exponents, edge_order = get_forest_exponents(n, roots_list)
        
        # Build Polyhedron
        # Use exact arithmetic
        verts = [vector(QQ, v) for v in exponents]
        P = Polyhedron(vertices=verts)
        
        # Get inequalities
        # form: A*x + b >= 0
        ieqs = P.inequality_generator()
        
        facets_desc = []
        for ieq in ieqs:
            # ieq is a vector (b, a1, ..., am)
            b = ieq[0]
            a = ieq[1:]
            
            # Check if it's a coordinate boundary: x_e >= 0 or x_e <= 1 (since binary)
            # x_e >= 0  => a has one 1, rest 0. b=0.
            # x_e <= 1  => -x_e >= -1 => a has one -1, rest 0. b=1.
            
            is_coord = False
            nz_indices = [i for i, val in enumerate(a) if val != 0]
            
            if len(nz_indices) == 1:
                idx = nz_indices[0]
                edge = edge_order[idx]
                val = a[idx]
                
                if val > 0 and b == 0:
                    facets_desc.append(f"z_{edge} -> 0")
                    if edge not in edge_coverage: edge_coverage[edge] = []
                    edge_coverage[edge].append(roots_list)
                    is_coord = True
                elif val < 0:
                    # x_e <= b/|val|
                    facets_desc.append(f"z_{edge} <= {b/abs(val)}")
                    is_coord = True
            
            if not is_coord:
                # Describe mixed facet
                # Just list edges involved
                terms = []
                for i in nz_indices:
                    terms.append(f"{a[i]}*x_{edge_order[i]}")
                facets_desc.append(f"Mixed: {' + '.join(terms)} + {b} >= 0")
                
        report_lines.append(f"## Roots {roots_list}")
        report_lines.append(f"- Vertices: {len(verts)}")
        report_lines.append(f"- Facets: {len(facets_desc)}")
        for f in facets_desc:
            report_lines.append(f"  - {f}")
        report_lines.append("")
        
    # Summary
    report_lines.append("## Edge Coordinate Boundary Coverage")
    report_lines.append("Which edges 'e' appear as simple facets x_e >= 0 (z_e -> 0)?")
    
    all_edges = []
    for i in range(n):
        for j in range(i+1, n):
            all_edges.append((i,j))
            
    for e in all_edges:
        if e in edge_coverage:
            roots_str = ", ".join([str(r) for r in edge_coverage[e]])
            report_lines.append(f"- z_{e}: Covered in {len(edge_coverage[e])} charts. (e.g. {roots_str})")
        else:
            report_lines.append(f"- z_{e}: **NOT COVERED** in any chart as a simple boundary.")
            
    # Write report
    with open("docs/root_patch_atlas.md", "w") as f:
        f.write("\n".join(report_lines))
        
    print("Atlas generated at docs/root_patch_atlas.md")

if __name__ == "__main__":
    analyze_root_sweep()




