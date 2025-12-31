import itertools
from sage.all import Graph

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

def get_forest_polynomial(forests, weights_func, C_func):
    """
    Computes sum_{F} prod_{(i,j)in F} (- w_ij C_i C_j)
    Note: The sign convention in expansion of det(L) usually involves -w for off-diagonals.
    And L_ii = sum w_ik ...
    Standard result: det(L^(R)) = sum_{F} prod_{(i,j) in F} (weight of edge)
    The weight of edge in L is w_ij C_i C_j.
    Wait, L_ij = - w_ij C_i C_j.
    Does the determinant expansion carry the minus signs?
    For Laplacian matrix where off-diagonals are negative, the determinant (All-Minors MTT)
    is exactly sum over forests of Product( - L_uv ).
    Since -L_uv = w_uv C_u C_v, the terms are positive.
    So we sum prod(w_ij C_i C_j).
    """
    total = 0
    for edges in forests:
        term = 1
        for u, v in edges:
            w = weights_func(u, v)
            c_u = C_func(u)
            c_v = C_func(v)
            term *= (w * c_u * c_v)
        total += term
    return total

