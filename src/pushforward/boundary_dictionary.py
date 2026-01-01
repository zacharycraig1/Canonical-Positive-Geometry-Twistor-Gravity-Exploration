"""
Boundary Dictionary for the Forest Polytope
Maps facets of the polytope to physical singularities (soft limits, collinear limits, factorization channels).
"""

class BoundaryDictionary:
    def __init__(self, n, roots, facets):
        """
        n: number of particles
        roots: list of root indices
        facets: list of facets (from Polytope.inequalities())
        """
        self.n = n
        self.roots = roots
        self.facets = facets
        self.mapping = {} # facet_index -> description
        
    def classify_facets(self):
        """
        Analyze each facet inequality and map it to a physical limit.
        Inequalities are of form: sum_{e in S} x_e <= rhs
        
        Types:
        1. x_e >= 0  (Lower bound) -> e -> 0?
        2. x_e <= 1  (Upper bound) -> e -> 1?
        3. Subsets constraints.
        """
        pass
        
    def get_physical_limit(self, facet_index):
        return self.mapping.get(facet_index, "Unknown")



