import sys
import os
from sage.all import *

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

from posgeom.forest_polytope import get_forest_exponents
from posgeom.canonical_polytope import eval_canonical_form_dual
from posgeom.physics_map import eval_edge_vars_from_spinors

def generate_random_kinematics(n):
    # Random spinors in Complex Field
    CC = ComplexField(200) # High precision
    lambdas = {i: vector(CC, [CC.random_element(), CC.random_element()]) for i in range(n)}
    tildes = {i: vector(CC, [CC.random_element(), CC.random_element()]) for i in range(n)}
    x = vector(CC, [CC.random_element(), CC.random_element()])
    y = vector(CC, [CC.random_element(), CC.random_element()])
    return lambdas, tildes, x, y

def run_test():
    print("Running W = (1, log|z|) diagnostic test for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Generate Kinematics
    lambdas, tildes, x, y = generate_random_kinematics(n)
    
    # 2. Compute z_ij
    # Need to adapt physics_map logic slightly if not directly importable or if it needs modification
    # But we can use the imported function if it works.
    # Note: physics_map uses simple dictionaries.
    
    # Re-implementing briefly to ensure control over types/precision
    z_map = {}
    CC = ComplexField(200)
    
    def bracket(l1, l2): return l1[0]*l2[1] - l1[1]*l2[0]
    
    # Weights
    C = {}
    for i in range(n):
        C[i] = bracket(lambdas[i], x) * bracket(lambdas[i], y)
        
    # Edges
    z_vals = {}
    
    # Get edge order from forest_polytope
    exponents, edge_order = get_forest_exponents(n, roots)
    
    print(f"Number of forest vertices: {len(exponents)}")
    print(f"Number of edges (dim): {len(edge_order)}")
    
    z_vector = []
    
    for (u, v) in edge_order:
        # z_{uv}
        ang = bracket(lambdas[u], lambdas[v])
        sq = bracket(tildes[u], tildes[v])
        val = (sq / ang) * C[u] * C[v]
        z_vals[(u, v)] = val
        z_vector.append(val)
        
    # 3. Compute W = (1, y) with y = log|z|
    # Note: canonical_form usually assumes W . Z = 1 + y . a_F
    # If we use log|z|, we are testing if the form behaves like the amplitude in "log variables".
    # This is heuristic.
    
    y_vector = [log(abs(val)) for val in z_vector]
    W = [1.0] + [float(val) for val in y_vector]
    
    # 4. Evaluate Canonical Form
    # Need to handle lower dimensional polytopes by projecting
    
    verts_qq = [vector(QQ, v) for v in exponents]
    P = Polyhedron(vertices=verts_qq)
    dim = P.dim()
    ambient_dim = len(exponents[0])
    print(f"Polyhedron dimension: {dim} (Ambient: {ambient_dim})")
    
    if dim < ambient_dim:
        print("Projecting to intrinsic dimension...")
        # 1. Pick origin
        v0 = verts_qq[0]
        
        # 2. Find basis of subspace
        # We can use the difference vectors
        diffs = [v - v0 for v in verts_qq[1:]]
        V_space = VectorSpace(QQ, ambient_dim)
        subspace = V_space.subspace(diffs)
        basis = subspace.basis()
        
        if len(basis) != dim:
             print(f"Warning: Basis length {len(basis)} != dim {dim}")
             
        # 3. Project vertices
        # v = v0 + sum c_i b_i
        # We need to find c_i for each vertex.
        # Construct matrix of basis vectors (columns)
        B_matrix = Matrix(basis).transpose()
        # Solve B * c = v - v0
        
        new_verts = []
        for v in verts_qq:
            diff = v - v0
            try:
                c = B_matrix.solve_right(diff)
                new_verts.append(list(c))
            except ValueError:
                print("Vertex not in subspace!")
                return
                
        # 4. Project W
        # W_old . (1, v) = w0 + w_vec . (v0 + B*c)
        #                = (w0 + w_vec . v0) + (w_vec . B) . c
        # So W_new = [w0 + w_vec . v0, (w_vec . B)_1, ...]
        
        w0 = W[0]
        w_vec = vector(RDF, W[1:])
        v0_rdf = vector(RDF, v0)
        B_rdf = Matrix(RDF, basis).transpose()
        
        w_new_0 = w0 + w_vec.dot_product(v0_rdf)
        w_new_rest = w_vec * B_rdf # Row vector * Matrix = Row vector
        
        W_new = [w_new_0] + list(w_new_rest)
        
        # Triangulate using exact arithmetic
        print("Triangulating polytope (exact)...")
        verts_qq_proj = [vector(QQ, v) for v in new_verts]
        # We need to manually triangulate because eval_canonical_form_dual does it on the input vertices
        # But we want to pass RDF vertices to eval, so we must pre-calculate triangulation on QQ vertices
        
        from posgeom.canonical_polytope import triangulate_polytope
        triangulation = triangulate_polytope(verts_qq_proj)
        print(f"Triangulation has {len(triangulation)} simplices")
        
        verts_final = [vector(RDF, v) for v in new_verts]
        W_final = vector(RDF, W_new)
        
    else:
        print("Triangulating polytope (exact)...")
        triangulation = triangulate_polytope(verts_qq)
        print(f"Triangulation has {len(triangulation)} simplices")
        
        verts_final = [vector(RDF, v) for v in exponents]
        W_final = vector(RDF, W)
    
    try:
        omega_val = eval_canonical_form_dual(W_final, verts_final, triangulation=triangulation)
        print(f"Omega(W) = {omega_val}")
    except Exception as e:
        print(f"Error evaluating canonical form: {e}")
        return

    # 5. Compute Target (Forest Polynomial / MHV Amplitude)
    # Target = Sum_F Prod_{e in F} z_e
    
    forest_poly_val = CC(0)
    for exp_vec in exponents:
        term = CC(1)
        for i, power in enumerate(exp_vec):
            if power == 1:
                term *= z_vector[i]
        forest_poly_val += term
        
    # Normalize amplitude
    # M_MHV ~ <xy>^8 * F(z) / (Prod C_k^2 ...)
    # We just compare Ratio = Omega / F(z) to see if it's constant-ish or wild.
    
    print(f"Forest Polynomial F(z) = {forest_poly_val}")
    print(f"|F(z)| = {abs(forest_poly_val)}")
    
    ratio = abs(omega_val) / abs(forest_poly_val)
    print(f"Ratio |Omega| / |F(z)| = {ratio}")
    
    # Also compare with naive W = (1, z)
    # This is expected to be garbage, but let's check.
    W_naive = [1.0] + [abs(z) for z in z_vector] # Using abs(z) to keep it real? Or should use complex?
    # Canonical form code usually assumes real polytope/dual.
    # If we use complex W, we need complex volume.
    # Let's stick to the log test which is the requested diagnostic.
    
    print("\nDiagnostic Conclusion:")
    print("If Ratio is constant across runs, there is a link.")
    print("If Ratio varies wildly, the identification is wrong (as expected).")

if __name__ == "__main__":
    run_test()

