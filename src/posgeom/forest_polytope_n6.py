from sage.all import *
import itertools

def get_forest_polynomial_n6(roots=[0,1,2]):
    """
    Optimized forest polynomial generator for n=6.
    Uses Prüfer sequences or direct enumeration to avoid large graph ops.
    
    Args:
        roots (list): List of 3 roots.
        
    Returns:
        Sage polynomial.
    """
    n = 6
    if len(roots) != 3:
        raise ValueError("Roots must be size 3")
        
    # Variables z_ij
    edge_vars = []
    for i in range(n):
        for j in range(i + 1, n):
            edge_vars.append(f"z_{i}_{j}")
    R = PolynomialRing(QQ, edge_vars)
    z = R.gens_dict()
    
    # Enumeration Strategy:
    # A 3-rooted forest on 6 vertices has 3 trees.
    # Vertices {0..5}. Roots {0,1,2}.
    # Non-roots {3,4,5}.
    #
    # Partitions of {3,4,5} into 3 sets (possibly empty) assigned to roots 0,1,2.
    #
    # Partition sizes (sum=3):
    # 1. 3, 0, 0 (3 perms of sizes)
    # 2. 2, 1, 0 (6 perms)
    # 3. 1, 1, 1 (1 perm)
    
    non_roots = [r for r in range(n) if r not in roots]
    
    poly = R(0)
    count = 0
    
    # Iterate over assignments of non-roots to roots
    # Each non-root chooses a root owner. 3^3 = 27 assignments.
    for assignment in itertools.product(roots, repeat=len(non_roots)):
        # Build components
        comps = {r: [r] for r in roots}
        for i, owner in enumerate(assignment):
            comps[owner].append(non_roots[i])
            
        # For each component, sum over spanning trees
        term_product = R(1)
        
        for r in roots:
            nodes = comps[r]
            if len(nodes) == 1:
                # Single node tree (root only) -> weight 1
                pass
            else:
                # Cayley's formula / Spanning trees of complete graph on 'nodes'
                # We need explicit edges to form monomials.
                # Generate spanning trees of K_{|nodes|}
                
                # Optimized for small sizes:
                k = len(nodes)
                # k=2: 1 tree (edge)
                # k=3: 3 trees (2 edges)
                # k=4: 16 trees (3 edges)
                
                if k == 2:
                    u, v = nodes[0], nodes[1]
                    if u > v: u, v = v, u
                    term_product *= z[f"z_{u}_{v}"]
                elif k == 3:
                    # Sum of 3 trees
                    u, v, w = nodes
                    # Edges: (uv, vw), (uv, uw), (vw, uw)
                    # Sort pairs
                    e1 = tuple(sorted((u,v)))
                    e2 = tuple(sorted((v,w)))
                    e3 = tuple(sorted((u,w)))
                    
                    z1 = z[f"z_{e1[0]}_{e1[1]}"]
                    z2 = z[f"z_{e2[0]}_{e2[1]}"]
                    z3 = z[f"z_{e3[0]}_{e3[1]}"]
                    
                    term_product *= (z1*z2 + z1*z3 + z2*z3)
                elif k == 4:
                    # 16 trees. 
                    # Use Sage Graph for robustness (k=4 is small enough)
                    sub_edges = list(itertools.combinations(nodes, 2))
                    G = Graph([nodes, sub_edges])
                    trees = G.spanning_trees()
                    
                    comp_sum = R(0)
                    for t in trees:
                        t_term = R(1)
                        for u, v, _ in t.edges():
                            if u > v: u, v = v, u
                            t_term *= z[f"z_{u}_{v}"]
                        comp_sum += t_term
                    term_product *= comp_sum
                    
        poly += term_product
        
    return poly







