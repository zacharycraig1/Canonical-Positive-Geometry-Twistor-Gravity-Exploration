import numpy as np
from sage.all import *

class CanonicalFormEvaluator:
    """
    Evaluates the canonical form of a polytope P on dual vector W.
    Omega_P(W) = sum_{triangulation} sgn(Delta) * Omega_Delta(W)
    
    where Omega_Delta(W) = det(Z_0...Z_d) / prod(W.Z_i)
    """
    
    @staticmethod
    def eval_simplex(vertices, W):
        """
        Evaluate canonical form for a simplex defined by vertices.
        vertices: list of vectors (or Matrix columns) [Z_0, ..., Z_d]
        W: dual vector (same dimension as Z_i)
        """
        dim = len(vertices) - 1
        # Check dimensions
        if len(W) != len(vertices[0]):
            raise ValueError(f"Dimension mismatch: W has {len(W)}, Z has {len(vertices[0])}")
            
        # Construct matrix Z = [Z_0, ..., Z_d]
        Z_mat = Matrix(vertices).transpose() # Columns are vertices
        
        # Numerator: det(Z)
        # Note: Sage Matrix determinant
        num = Z_mat.det()
        
        # Denominator: prod(W . Z_i)
        denom = 1
        for Z_i in vertices:
            dot_prod = sum(w*z for w, z in zip(W, Z_i))
            if dot_prod == 0:
                raise ValueError("Pole encountered: W . Z_i = 0")
            denom *= dot_prod
            
        return num / denom

    @staticmethod
    def eval_polytope(P, W, method='triangulation'):
        """
        Evaluate canonical form for a polytope P.
        P: Sage Polyhedron object
        W: dual vector
        method: 'triangulation' (default)
        """
        # 1. Handle Dimension Deficiency (Projection)
        # Check if the ambient dimension matches the intrinsic dimension
        # P.vertices() are in R^N. P.dim() is d.
        # If N > d, we must project to R^d to define the volume form consistently.
        
        vertices = list(P.vertices())
        vecs = [v.vector() for v in vertices]
        N = len(vecs[0])
        d = P.dim()
        
        # W has dimension N+1 (homogeneous dual)
        if len(W) != N + 1:
            raise ValueError(f"W must have dimension {N+1} (ambient dim + 1), got {len(W)}")

        if d < N:
            # PROJECT TO AFFINE SUBSPACE
            # 1. Pick origin v0
            v0 = vecs[0]
            
            # 2. Find basis for subspace V = span(v - v0)
            shifted_vecs = [v - v0 for v in vecs]
            # Use Sage's VectorSpace to find basis
            V_space = VectorSpace(QQ, N)
            subspace = V_space.subspace(shifted_vecs)
            basis = subspace.basis() # List of d vectors B_1, ..., B_d
            
            if len(basis) != d:
                # This check ensures numerical stability / correct dimension
                raise ValueError(f"Basis size {len(basis)} does not match P.dim() {d}")
                
            # 3. Map vertices to coordinates c in R^d
            # v = v0 + sum c_i B_i
            # We need to solve for c_i.
            # Construct matrix B with columns B_i. (N x d)
            # (v - v0) = B * c
            # c = (B^T B)^-1 B^T (v - v0)  (or use solve_right)
            
            B_mat = Matrix(basis).transpose() # N x d
            
            projected_vecs = []
            for v in vecs:
                diff = v - v0
                # solve B_mat * c = diff
                try:
                    c = B_mat.solve_right(diff)
                    projected_vecs.append(c)
                except ValueError:
                    raise ValueError("Vertex not in affine span!")
            
            # 4. Map W to W_proj in (R^d)*
            # W_proj = [ W.(1, v0), W.(0, B_1), ..., W.(0, B_d) ]
            # W is (N+1)-dim: (w_0, w_1...w_N)
            # Dual pairing: W . (1, v) = w_0 + sum w_k v_k
            
            # Component 0: W evaluated at origin v0
            # (1, v0)
            Z0 = vector(QQ, [1] + list(v0))
            w_proj_0 = W.dot_product(Z0)
            
            # Component i: W evaluated at direction B_i (homogeneous part 0)
            w_proj_rest = []
            for b in basis:
                # (0, b)
                Zb = vector(QQ, [0] + list(b))
                val = W.dot_product(Zb)
                w_proj_rest.append(val)
                
            W_proj = vector(QQ, [w_proj_0] + w_proj_rest)
            
            # Update variables for evaluation
            # Vertices are now in R^d
            vecs = projected_vecs
            # W is in R^{d+1}
            W = W_proj
            
            # Update N to be d
            N = d
            
            # Note: We need to reconstruct P for triangulation?
            # Yes, we need to triangulate the projected polytope.
            P_proj = Polyhedron(vertices=vecs)
            P = P_proj

        if method == 'triangulation':
            # Use Sage's triangulation
            try:
                triangulation = P.triangulation()
            except AttributeError:
                # Fallback using PointConfiguration for backends like PPL
                from sage.geometry.triangulation.point_configuration import PointConfiguration
                pc = PointConfiguration(P.vertices())
                triangulation = pc.triangulate()

            # Homogenize vertices for evaluation: Z = (1, v)
            Z_vecs = [vector(QQ, [1] + list(v)) for v in vecs]

            # Determine a reference W_ref in the dual cone
            dim = len(vecs[0]) # Should be d
            W_ref = vector(QQ, [1]*(dim+1))
            
            # Perturb W_ref if invalid
            if any(W_ref.dot_product(z) == 0 for z in Z_vecs):
                W_ref = vector(QQ, [1]*(dim+1)) + vector(QQ, [QQ(1)/(i+10) for i in range(dim+1)])

            total_val = 0
            # P.triangulation() yields tuples of indices
            for simplex_indices in triangulation:
                simplex_verts = [Z_vecs[i] for i in simplex_indices]
                
                # Evaluate at W_ref to determine sign
                try:
                    ref_val = CanonicalFormEvaluator.eval_simplex(simplex_verts, W_ref)
                except ValueError:
                    # Retry with random W_ref
                    W_ref_rand = CanonicalFormEvaluator.get_valid_W(P, seed=42)
                    ref_val = CanonicalFormEvaluator.eval_simplex(simplex_verts, W_ref_rand)
                
                if ref_val == 0: continue 

                sign = 1 if ref_val > 0 else -1
                
                # Evaluate at actual W
                term = CanonicalFormEvaluator.eval_simplex(simplex_verts, W)
                total_val += sign * term
                
            return total_val
        else:
            raise NotImplementedError(f"Method {method} not implemented")

    @staticmethod
    def get_valid_W(P, seed=None):
        """
        Generate a random W that avoids poles (W . Z_i != 0).
        P should be the projected polytope if applicable.
        """
        if seed is not None:
            set_random_seed(seed)
            
        dim = P.dim() + 1 # Ambient dimension (projective)
        vertices = [vector(QQ, [1] + list(v)) for v in P.vertices()]
        
        while True:
            # Random integer vector
            W = vector(QQ, [randint(-10, 10) for _ in range(dim)])
            if all(W.dot_product(v) != 0 for v in vertices):
                return W
