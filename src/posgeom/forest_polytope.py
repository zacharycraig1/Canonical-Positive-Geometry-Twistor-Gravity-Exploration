import itertools
from sage.all import PolynomialRing, QQ, Graph

def enumerate_rooted_forests(n, roots):
    """
    Enumerates all spanning forests of K_n with |roots| components,
    where each component contains exactly one root from 'roots'.
    
    Returns a list of edge lists. Each edge is (u, v) with u < v.
    """
    non_roots = [i for i in range(n) if i not in roots]
    k = len(roots)
    
    forests = []
    
    # Iterate over all assignments of non_roots to roots
    # Each non_root i is assigned to a root r (index 0..k-1)
    for assignment in itertools.product(range(k), repeat=len(non_roots)):
        # Build components
        components = [ [roots[i]] for i in range(k) ]
        for idx, owner in enumerate(assignment):
            components[owner].append(non_roots[idx])
            
        # For each component, generate all spanning trees
        # Since it's a complete graph K_n restricted to the component vertices
        component_trees = []
        possible = True
        
        for comp in components:
            if len(comp) == 1:
                component_trees.append([[]]) # Single vertex tree has no edges
            else:
                # Generate spanning trees of K_{|comp|}
                # We can use Sage's Graph or Cayley's
                # For small sets, simpler to just use generic spanning tree generator
                subgraph_edges = list(itertools.combinations(comp, 2))
                G = Graph([comp, subgraph_edges])
                trees = G.spanning_trees()
                # Extract edges
                edges_list = []
                for t in trees:
                    edges_list.append(t.edges(labels=False, sort=True))
                component_trees.append(edges_list)
        
        # Cartesian product of trees from each component
        for forest_tuple in itertools.product(*component_trees):
            # Flatten edges
            full_forest_edges = []
            for tree_edges in forest_tuple:
                full_forest_edges.extend(tree_edges)
            forests.append(full_forest_edges)
            
    return forests

def get_forest_polynomial(n, roots):
    """
    Computes the forest polynomial F_{n,R}(z) for the complete graph K_n
    with respect to roots R.
    
    The polynomial is a sum over all spanning forests F where each component
    contains exactly one root from 'roots'.
    
    Monomial for forest F: prod_{(i,j) in E(F)} z_{ij}
    
    Args:
        n (int): Number of vertices (0 to n-1)
        roots (list): List of root vertices, e.g. [0, 1, 2]
        
    Returns:
        Sage polynomial in variables z_{ij} (for i<j)
    """
    # 1. Define variables z_{ij} for 0 <= i < j < n
    edge_vars = []
    for i in range(n):
        for j in range(i + 1, n):
            edge_vars.append(f"z_{i}_{j}")
            
    R = PolynomialRing(QQ, edge_vars)
    z = R.gens_dict()
    
    # 2. Enumerate rooted forests
    forests_edges = enumerate_rooted_forests(n, roots)
    
    # 3. Build Polynomial
    poly = R(0)
    for edges in forests_edges:
        term = R(1)
        for u, v in edges:
            if u > v: u, v = v, u
            term *= z[f"z_{u}_{v}"]
        poly += term
        
    return poly

def get_forest_exponents(n, roots):
    """
    Returns the list of exponent vectors for the forest polynomial.
    Each vector corresponds to a forest.
    Ordering of edges is canonical: (0,1), (0,2), ..., (n-2, n-1).
    """
    # Canonical edge order
    edges_ordered = []
    for i in range(n):
        for j in range(i + 1, n):
            edges_ordered.append((i, j))
            
    forests_edges = enumerate_rooted_forests(n, roots)
    
    exponents = []
    dim = len(edges_ordered)
    edge_to_idx = {e: k for k, e in enumerate(edges_ordered)}
    
    for forest in forests_edges:
        vec = [0] * dim
        for u, v in forest:
            if u > v: u, v = v, u
            idx = edge_to_idx[(u, v)]
            vec[idx] = 1
        exponents.append(vec)
        
    return exponents, edges_ordered
