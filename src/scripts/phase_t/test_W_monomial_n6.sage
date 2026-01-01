import sys
import os
from sage.all import *

sys.path.append(os.path.join(os.getcwd(), 'src'))

from posgeom.forest_polytope import get_forest_exponents
from posgeom.canonical_polytope import eval_canonical_form_dual, triangulate_polytope

def generate_random_kinematics(n):
    CC = ComplexField(200)
    lambdas = {i: vector(CC, [CC.random_element(), CC.random_element()]) for i in range(n)}
    tildes = {i: vector(CC, [CC.random_element(), CC.random_element()]) for i in range(n)}
    x = vector(CC, [CC.random_element(), CC.random_element()])
    y = vector(CC, [CC.random_element(), CC.random_element()])
    return lambdas, tildes, x, y

def run_test():
    print("Running W = (z^a_F) Monomial Map test for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Kinematics
    lambdas, tildes, x, y = generate_random_kinematics(n)
    
    def bracket(l1, l2): return l1[0]*l2[1] - l1[1]*l2[0]
    C = {}
    for i in range(n): C[i] = bracket(lambdas[i], x) * bracket(lambdas[i], y)
    
    # 2. Edges and z_vals
    exponents, edge_order = get_forest_exponents(n, roots)
    z_map = {}
    z_vector = [] # indexed by edge_order
    
    for (u, v) in edge_order:
        ang = bracket(lambdas[u], lambdas[v])
        sq = bracket(tildes[u], tildes[v])
        val = (sq / ang) * C[u] * C[v]
        z_vector.append(val)
        
    # 3. Construct Monomial W
    # W_F = product z_e for e in F
    W_vals = []
    for exp in exponents:
        term = 1
        for i, p in enumerate(exp):
            if p == 1:
                term *= z_vector[i]
        W_vals.append(term)
        
    # We need real part for canonical form eval usually?
    # Or we can try abs(W) test first to avoid complex triangulation issues?
    # Let's try W = abs(W_vals) first.
    W_abs = [abs(w) for w in W_vals]
    
    # 4. Canonical Form
    # Project polytope
    verts_qq = [vector(QQ, v) for v in exponents]
    v0 = verts_qq[0]
    diffs = [v - v0 for v in verts_qq[1:]]
    V_space = VectorSpace(QQ, len(exponents[0]))
    subspace = V_space.subspace(diffs)
    basis = subspace.basis()
    
    # Project W
    # W . Z = sum W_F * 1 (homogenized is just 1... wait)
    # The homogenized coord Z_F = (1, a_F).
    # Dual W has dim d+1.
    # But here W has dim N (number of vertices).
    # WAIT.
    # Omega(W) takes W in dual of EMBEDDING space (dim d+1).
    # The W I constructed has size 108 (number of vertices).
    # This is NOT a dual vector.
    # This is a set of coefficients defining the polynomial $F(u) = \sum W_F u^F$.
    # But `eval_canonical_form_dual` expects a linear form $L(Z) = W \cdot Z$.
    # It evaluates $\Omega(L)$.
    # This means it assumes the polynomial is linear: $L(Z)$.
    # BUT our polynomial $F(u)$ is NOT linear in $u$. It is linear in $X = u^F$.
    # The canonical form of the polytope $P$ is a form on the dual of the space containing $P$.
    # $P \subset \mathbb{R}^{15}$.
    # Dual is $\mathbb{R}^{16}$.
    # But my "Monomial W" hypothesis gives 108 numbers.
    # This mismatch confirms: The Amplitude is NOT $\Omega_P(W)$.
    # The Amplitude is related to the coefficients $W_F$.
    # But $\Omega_P$ depends on a vector of size 16.
    
    print("MISMATCH DETECTED:")
    print(f"W vector from monomials has size {len(W_vals)} (Number of Forests).")
    print(f"Polytope dual space has dimension {len(exponents[0]) + 1} (16).")
    print("Cannot evaluate Omega_P at Monomial W directly.")
    print("Conclusion: The map is not simply W = z.")
    
    # Check if W can be represented as a linear form?
    # i.e. is there a L such that L(a_F) = W_F?
    # This would mean log(z) relation?
    # If W_F = exp(L . a_F) then log W_F = L . a_F.
    # log(prod z_e) = sum log z_e.
    # So L should be (0, log z_1, ..., log z_m).
    # This returns us to the Log Z test.
    # And we already tested Log Z and it failed.
    
    print("Re-evaluating Log Z test result interpretation...")
    # If log Z failed, then the relation is not L . a_F.
    # The amplitude is not the canonical form evaluated at L = log z.
    
    print("Test aborted due to dimensional mismatch.")

if __name__ == "__main__":
    run_test()



