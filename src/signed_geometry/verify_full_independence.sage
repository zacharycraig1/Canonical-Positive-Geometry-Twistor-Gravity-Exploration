#!/usr/bin/env sage
"""
Verify Full Reference Independence of Gravity Amplitude

The forest sum det(L~) depends on reference spinors.
But the FULL AMPLITUDE (with proper normalization) must be independent.

Formula:
    M = (-1)^(n-1) × ⟨ab⟩^8 × det(L~) / (prod_k C_k^2 × N_R)
    
The C_k factors in numerator and denominator should cancel.
"""
from sage.all import *
import itertools
import sys
import os

sys.path.insert(0, os.getcwd())

load('src/spinor_sampling.sage')


def ang_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def sq_bracket(tilde_lambdas, i, j):
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def enumerate_rooted_forests(n, roots):
    """Enumerate all spanning forests with specified roots."""
    num_roots = len(roots)
    num_edges = n - num_roots
    roots_set = set(roots)
    
    all_edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    
    for edges in itertools.combinations(all_edges, num_edges):
        adj = {i: [] for i in range(n)}
        for u, v in edges:
            adj[u].append(v)
            adj[v].append(u)
        
        visited = set()
        components = []
        valid = True
        
        for i in range(n):
            if i not in visited:
                stack = [i]
                visited.add(i)
                comp = [i]
                root_count = 1 if i in roots_set else 0
                
                while stack:
                    curr = stack.pop()
                    for neighbor in adj[curr]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)
                            comp.append(neighbor)
                            if neighbor in roots_set:
                                root_count += 1
                
                components.append(comp)
                if root_count != 1:
                    valid = False
                    break
        
        if valid and len(components) == num_roots:
            yield edges


def compute_full_amplitude(lambdas, tilde_lambdas, x_spinor, y_spinor, roots=(0, 1, 2)):
    """
    Compute the fully normalized gravity amplitude.
    
    This should be independent of reference spinors.
    """
    n = len(lambdas)
    non_roots = [i for i in range(n) if i not in roots]
    
    def ang_with_ref(lam, ref):
        return lam[0] * ref[1] - lam[1] * ref[0]
    
    # Compute C_i = ⟨i,x⟩⟨i,y⟩
    C = {i: ang_with_ref(lambdas[i], x_spinor) * ang_with_ref(lambdas[i], y_spinor) 
         for i in range(n)}
    
    # Compute w_ij = [ij]/⟨ij⟩
    w = {}
    for i in range(n):
        for j in range(i+1, n):
            ang = ang_bracket(lambdas, i, j)
            sq = sq_bracket(tilde_lambdas, i, j)
            if ang != 0:
                w[(i, j)] = sq / ang
    
    # Compute det(L~) = sum over forests
    det_L = QQ(0)
    for forest in enumerate_rooted_forests(n, roots):
        edges = list(forest)
        weight = QQ(1)
        for u, v in edges:
            if (u, v) not in w:
                continue
            weight *= (-w[(u, v)] * C[u] * C[v])
        det_L += weight
    
    # Normalization: product of C_k^2 for non-roots
    C_denom = QQ(1)
    for k in non_roots:
        C_denom *= C[k]**2
    
    # Parke-Taylor normalization for roots
    pt_norm = QQ(1)
    for i in range(len(roots)):
        j = (i + 1) % len(roots)
        pt_norm *= ang_bracket(lambdas, roots[i], roots[j])
    pt_norm = pt_norm**2
    
    # Helicity factor
    ang_01 = ang_bracket(lambdas, 0, 1)
    
    # Full amplitude
    # M = (-1)^(n-1) × ⟨ab⟩^8 × det(L~) / (C_denom × pt_norm)
    if C_denom == 0 or pt_norm == 0:
        return None, None, None
    
    M = (-1)**(n-1) * ang_01**8 * det_L / (C_denom * pt_norm)
    
    return M, det_L, C_denom


def verify_independence():
    """Verify that the full amplitude is reference-independent."""
    print("=" * 70)
    print("FULL AMPLITUDE REFERENCE INDEPENDENCE")
    print("=" * 70)
    
    n = 6
    roots = (0, 1, 2)
    
    for sample in range(5):
        print(f"\n--- Sample {sample} ---")
        
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 97)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        # Try different reference spinors
        reference_pairs = [
            (vector(QQ, [1, 2]), vector(QQ, [3, 1])),
            (vector(QQ, [1, 0]), vector(QQ, [0, 1])),
            (vector(QQ, [5, 7]), vector(QQ, [11, 3])),
            (vector(QQ, [-2, 5]), vector(QQ, [7, -3])),
        ]
        
        amplitudes = []
        for x_spinor, y_spinor in reference_pairs:
            M, det_L, C_denom = compute_full_amplitude(
                lambdas, tilde_lambdas, x_spinor, y_spinor, roots
            )
            
            if M is not None:
                amplitudes.append(M)
                print(f"  Ref {list(x_spinor)}: M = {float(M):.6e}, det_L = {float(det_L):.2e}, C_denom = {float(C_denom):.2e}")
        
        # Check if all amplitudes are equal
        if len(amplitudes) >= 2:
            base = amplitudes[0]
            all_equal = all(abs(float(M - base)) / abs(float(base)) < 1e-10 
                           if base != 0 else M == 0 
                           for M in amplitudes)
            
            if all_equal:
                print("  ✓ All amplitudes EQUAL (reference independent)")
            else:
                ratios = [float(M / base) if base != 0 else 0 for M in amplitudes]
                print(f"  ✗ Amplitudes vary: ratios = {ratios}")


def compare_with_hodges():
    """Compare the forest formula with the Hodges determinant."""
    print("\n" + "=" * 70)
    print("COMPARISON WITH HODGES DETERMINANT")
    print("=" * 70)
    
    load('src/hodges.sage')
    load('src/kinematics_map.sage')
    
    n = 6
    roots = (0, 1, 2)
    
    for sample in range(5):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 83)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        # Compute forest amplitude
        M_forest, det_L, C_denom = compute_full_amplitude(
            lambdas, tilde_lambdas, x_spinor, y_spinor, roots
        )
        
        if M_forest is None:
            continue
        
        # Compute Hodges amplitude
        Z_list = []
        x_current = matrix(QQ, 2, 2, 0)
        for i in range(n):
            lam = lambdas[i]
            til = tilde_lambdas[i]
            mu_0 = x_current[0,0]*lam[0] + x_current[0,1]*lam[1]
            mu_1 = x_current[1,0]*lam[0] + x_current[1,1]*lam[1]
            Z_list.append(vector(QQ, [lam[0], lam[1], mu_0, mu_1]))
            p_matrix = matrix(QQ, 2, 2, [lam[0]*til[0], lam[0]*til[1], 
                                         lam[1]*til[0], lam[1]*til[1]])
            x_current += p_matrix
        
        twistor = MomentumTwistor(n=n, Z=Z_list, check_domain=False)
        twistor._compute_brackets()
        
        M_hodges, reason = hodges_6pt_mhv(twistor)
        
        if M_hodges is None:
            print(f"Sample {sample}: Hodges failed - {reason}")
            continue
        
        # Compare
        if M_hodges != 0:
            ratio = M_forest / M_hodges
            print(f"Sample {sample}: M_forest/M_hodges = {float(ratio):.6f}")


if __name__ == "__main__":
    verify_independence()
    compare_with_hodges()

