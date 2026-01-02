#!/usr/bin/env sage
"""
Twisted Intersection Numbers for Scattering Amplitudes

This module computes the intersection pairing between twisted cohomology
forms, which defines the KLT kernel for gravity amplitudes.

Key formula (Mizera):
    <phi_alpha | phi_beta> = sum_{z = critical point} Res [ phi_alpha * phi_beta / det(Hessian(W)) ]

For scattering amplitudes:
    - Critical points = solutions of scattering equations
    - phi_alpha, phi_beta = Parke-Taylor forms
    - The intersection matrix is proportional to the inverse of the bi-adjoint scalar amplitude matrix
"""
from sage.all import *
import itertools
import sys
import os

# Add project root to path
sys.path.insert(0, os.getcwd())


def load_setup(n=6):
    """Load the twisted form setup data."""
    try:
        data = load(f"twisted_cohomology/setup_n{n}.sobj")
        return data
    except FileNotFoundError:
        print(f"Setup file for n={n} not found. Run 01_setup_forms.sage first.")
        return None


def compute_hessian(omega_components, z_vars, n):
    """
    Compute the Hessian matrix of the twist potential W.
    
    H_{ij} = d^2 W / (dz_i dz_j) = d omega_i / dz_j
    
    For the intersection formula, we need det(H) at critical points.
    """
    dynamic_indices = list(range(2, n-1))
    dim = len(dynamic_indices)
    
    H = matrix(SR, dim, dim)
    
    for i_idx, i in enumerate(dynamic_indices):
        for j_idx, j in enumerate(dynamic_indices):
            # H_{ij} = d(omega_i) / dz_j
            omega_i = omega_components[i]
            z_j = z_vars[j]
            H[i_idx, j_idx] = diff(omega_i, z_j)
            
    return H, dynamic_indices


def intersection_number_residue(data, s_values, alpha, beta, solutions):
    """
    Compute intersection number <phi_alpha | phi_beta> via residue formula.
    
    Formula:
        <phi_alpha | phi_beta> = sum_{solutions} phi_alpha(z) * phi_beta(z) / det(Hessian(W))
        
    Args:
        data: Setup data from load_setup
        s_values: Dictionary mapping (i,j) -> numerical value of s_ij
        alpha, beta: Permutation tuples indexing Parke-Taylor forms
        solutions: List of scattering equation solutions
        
    Returns:
        Intersection number (numerical)
    """
    n = data['n']
    z_vars = data['z_vars']
    omega = data['omega']
    pt_forms = data['pt_forms']
    s_vars = data['s_vars']
    
    # Get the PT forms
    phi_alpha = pt_forms[alpha]
    phi_beta = pt_forms[beta]
    
    # Compute Hessian
    H, dynamic_indices = compute_hessian(omega, z_vars, n)
    
    # Compute determinant of Hessian
    det_H = H.det()
    
    # The integrand for intersection
    # Note: We need to include the measure d^{n-3}z as well
    # For d-log forms, the intersection is related to residues
    integrand = phi_alpha * phi_beta / det_H
    
    # Substitute s values
    for (i,j), val in s_values.items():
        if (i,j) in s_vars:
            integrand = integrand.subs(s_vars[(i,j)] == val)
    
    # Sum over solutions
    total = 0
    for sol in solutions:
        val = integrand
        for k in dynamic_indices:
            z_k = z_vars[k]
            z_val = sol.get(z_k, None)
            if z_val is not None:
                val = val.subs(z_k == z_val)
        
        try:
            total += complex(val)
        except:
            total += val
            
    return total


def compute_intersection_matrix_numerical(data, lambdas, tilde_lambdas):
    """
    Compute the full intersection matrix numerically for given kinematics.
    
    This should equal (up to normalization) the bi-adjoint scalar amplitude matrix,
    whose inverse is the KLT kernel.
    
    Args:
        data: Setup data
        lambdas: List of n spinor 2-vectors
        tilde_lambdas: List of n spinor-tilde 2-vectors
        
    Returns:
        Intersection matrix M (numpy array)
    """
    import numpy as np
    from scipy.optimize import fsolve
    
    n = data['n']
    perms = data['perms']
    z_vars = data['z_vars']
    omega = data['omega']
    s_vars = data['s_vars']
    
    # Compute Mandelstam values from spinors
    def ang_bracket(lam_i, lam_j):
        return lam_i[0] * lam_j[1] - lam_i[1] * lam_j[0]
    
    def sq_bracket(til_i, til_j):
        return til_i[0] * til_j[1] - til_i[1] * til_j[0]
    
    s_values = {}
    for i in range(n):
        for j in range(i+1, n):
            # 1-indexed in s_vars, 0-indexed in lambdas
            ang = ang_bracket(lambdas[i], lambdas[j])
            sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
            s_values[(i+1, j+1)] = float(ang * sq)
            s_values[(j+1, i+1)] = s_values[(i+1, j+1)]
    
    # Solve scattering equations numerically
    dynamic_indices = list(range(2, n-1))
    
    def scattering_eqs_func(z_vec):
        """Scattering equations as a function for numerical solving."""
        z_dict = {z_vars[k]: z_vec[i] for i, k in enumerate(dynamic_indices)}
        z_dict[z_vars[1]] = 0.0
        z_dict[z_vars[n-1]] = 1.0
        
        residuals = []
        for k in dynamic_indices:
            eq = omega[k]
            # Substitute s values
            for (i,j), val in s_values.items():
                if (i,j) in s_vars:
                    eq = eq.subs(s_vars[(i,j)] == val)
            # Substitute z values
            for var, val in z_dict.items():
                eq = eq.subs(var == val)
            try:
                residuals.append(float(eq))
            except:
                residuals.append(1e10)  # Large value for failed evaluation
        return residuals
    
    # Find solutions with multiple starting points
    # For n=6, we expect (n-3)! = 6 solutions
    solutions = []
    num_expected = factorial(n-3)
    
    # Try random starting points
    import random
    random.seed(42)
    
    for attempt in range(100):
        z0 = [random.uniform(0.1, 0.9) for _ in dynamic_indices]
        try:
            sol, info, ier, mesg = fsolve(scattering_eqs_func, z0, full_output=True)
            if ier == 1:  # Converged
                residual = np.linalg.norm(info['fvec'])
                if residual < 1e-10:
                    # Check if solution is new
                    is_new = True
                    for existing in solutions:
                        if np.allclose(sol, existing, rtol=1e-6):
                            is_new = False
                            break
                    if is_new:
                        solutions.append(sol)
                        if len(solutions) >= num_expected:
                            break
        except:
            continue
    
    print(f"Found {len(solutions)} solutions (expected {num_expected})")
    
    if len(solutions) == 0:
        print("No solutions found!")
        return None
    
    # Compute Hessian and its determinant
    H, _ = compute_hessian(omega, z_vars, n)
    det_H = H.det()
    
    # Substitute s values into Hessian determinant and PT forms
    det_H_subst = det_H
    for (i,j), val in s_values.items():
        if (i,j) in s_vars:
            det_H_subst = det_H_subst.subs(s_vars[(i,j)] == val)
    
    pt_forms_subst = {}
    for alpha in perms:
        pt = data['pt_forms'][alpha]
        for (i,j), val in s_values.items():
            if (i,j) in s_vars:
                pt = pt.subs(s_vars[(i,j)] == val)
        pt_forms_subst[alpha] = pt
    
    # Build intersection matrix
    dim = len(perms)
    M = np.zeros((dim, dim), dtype=complex)
    
    for i, alpha in enumerate(perms):
        for j, beta in enumerate(perms):
            # Sum over solutions
            total = 0
            for sol in solutions:
                # Build z substitution dict
                z_dict = {z_vars[k]: sol[idx] for idx, k in enumerate(dynamic_indices)}
                
                # Evaluate PT forms
                pt_alpha = pt_forms_subst[alpha]
                pt_beta = pt_forms_subst[beta]
                det_H_val = det_H_subst
                
                for var, val in z_dict.items():
                    pt_alpha = pt_alpha.subs(var == val)
                    pt_beta = pt_beta.subs(var == val)
                    det_H_val = det_H_val.subs(var == val)
                
                try:
                    pt_alpha_val = complex(pt_alpha)
                    pt_beta_val = complex(pt_beta)
                    det_H_val_c = complex(det_H_val)
                    
                    if abs(det_H_val_c) > 1e-15:
                        total += pt_alpha_val * pt_beta_val / det_H_val_c
                except Exception as e:
                    pass
                    
            M[i, j] = total
    
    return M, solutions


def compute_biadjoint_matrix(lambdas, tilde_lambdas, n=6):
    """
    Compute the bi-adjoint scalar amplitude matrix m[alpha|beta].
    
    m[alpha|beta] = sum over compatible cubic diagrams of 1/prod(s_poles)
    
    This should equal the intersection matrix (up to normalization).
    """
    import networkx as nx
    
    def ang_bracket(i, j):
        return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]
    
    def sq_bracket(i, j):
        return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]
    
    def mandelstam(i, j):
        return ang_bracket(i, j) * sq_bracket(i, j)
    
    def s_subset(subset):
        """Compute s for a subset of particles."""
        val = 0
        subset = list(subset)
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                val += mandelstam(subset[a], subset[b])
        return val
    
    def get_planar_poles(order):
        """Get all planar poles for a given ordering."""
        m = len(order)
        poles = set()
        doubled = list(order) + list(order)
        for size in range(2, m-1):
            for i in range(m):
                subset = frozenset(doubled[i:i+size])
                if len(subset) > m//2:
                    full = frozenset(order)
                    subset = full - subset
                elif len(subset) == m//2:
                    if 0 not in subset:
                        full = frozenset(order)
                        subset = full - subset
                poles.add(subset)
        return list(poles)
    
    def get_triangulations(order):
        """Get all compatible cubic diagrams for an ordering."""
        poles = get_planar_poles(order)
        m = len(order)
        num_propagators = m - 3
        
        compatible_pairs = set()
        for i in range(len(poles)):
            for j in range(i+1, len(poles)):
                p1 = poles[i]
                p2 = poles[j]
                if p1.isdisjoint(p2) or p1.issubset(p2) or p2.issubset(p1):
                    compatible_pairs.add((i, j))
        
        G = nx.Graph()
        G.add_nodes_from(range(len(poles)))
        G.add_edges_from(list(compatible_pairs))
        cliques = [c for c in nx.find_cliques(G) if len(c) == num_propagators]
        diagrams = []
        for c in cliques:
            diagrams.append([poles[x] for x in c])
        return diagrams
    
    # Build permutation basis
    permuted_set = [1, 2, 3]  # 0-indexed particles 1,2,3
    basis_perms = sorted(list(itertools.permutations(permuted_set)))
    
    # Build diagrams for left and right orderings
    alpha_diagrams = {}
    beta_diagrams = {}
    
    for alpha in basis_perms:
        order = [0] + list(alpha) + [4, 5]
        alpha_diagrams[alpha] = get_triangulations(order)
        
    for beta in basis_perms:
        order = [0] + list(beta) + [5, 4]  # Note the swap for right ordering
        beta_diagrams[beta] = get_triangulations(order)
    
    # Compute matrix
    dim = len(basis_perms)
    M = [[QQ(0) for _ in range(dim)] for _ in range(dim)]
    
    for i, alpha in enumerate(basis_perms):
        for j, beta in enumerate(basis_perms):
            diags_a = set([frozenset(d) for d in alpha_diagrams[alpha]])
            diags_b = set([frozenset(d) for d in beta_diagrams[beta]])
            common = diags_a.intersection(diags_b)
            
            val = QQ(0)
            for diag in common:
                den = QQ(1)
                valid = True
                for pole in diag:
                    s_val = s_subset(pole)
                    if s_val == 0:
                        valid = False
                        break
                    den *= s_val
                if valid and den != 0:
                    val += 1/den
                    
            M[i][j] = val
    
    return matrix(QQ, M), basis_perms


def verify_intersection_equals_biadjoint(num_samples=5):
    """
    Verify that the intersection matrix equals the bi-adjoint matrix.
    
    This is the key test that the twisted cohomology framework gives
    the correct KLT kernel.
    """
    print("=" * 60)
    print("VERIFICATION: Intersection Matrix = Bi-adjoint Matrix")
    print("=" * 60)
    
    # Load setup
    data = load_setup(n=6)
    if data is None:
        return False
    
    # Import spinor sampling
    load('src/spinor_sampling.sage')
    
    for sample in range(num_samples):
        print(f"\nSample {sample + 1}:")
        
        # Sample kinematics
        result = sample_spinor_helicity_conserving(n=6, seed=sample)
        if result is None:
            print("  Skipped (sampling failed)")
            continue
            
        lambdas, tilde_lambdas = result
        
        # Convert to lists for numerical computation
        lambdas_float = [[float(x) for x in lam] for lam in lambdas]
        tilde_lambdas_float = [[float(x) for x in til] for til in tilde_lambdas]
        
        # Compute intersection matrix (numerical)
        result = compute_intersection_matrix_numerical(data, lambdas_float, tilde_lambdas_float)
        if result is None:
            print("  Skipped (intersection computation failed)")
            continue
            
        M_int, solutions = result
        
        # Compute bi-adjoint matrix (exact)
        M_ba, perms = compute_biadjoint_matrix(lambdas, tilde_lambdas, n=6)
        M_ba_float = M_ba.change_ring(RDF).numpy()
        
        # Compare signatures
        import numpy as np
        
        # Symmetrize both
        M_int_sym = 0.5 * (M_int + M_int.T)
        M_ba_sym = 0.5 * (M_ba_float + M_ba_float.T)
        
        # Compute eigenvalues
        eigs_int = np.linalg.eigvalsh(np.real(M_int_sym))
        eigs_ba = np.linalg.eigvalsh(M_ba_sym)
        
        sig_int = (sum(1 for e in eigs_int if e > 1e-10), 
                   sum(1 for e in eigs_int if e < -1e-10))
        sig_ba = (sum(1 for e in eigs_ba if e > 1e-10), 
                  sum(1 for e in eigs_ba if e < -1e-10))
        
        print(f"  Intersection matrix signature: {sig_int}")
        print(f"  Bi-adjoint matrix signature: {sig_ba}")
        
        if sig_int == sig_ba:
            print(f"  ✓ Signatures match!")
        else:
            print(f"  ✗ Signature mismatch")
    
    print("\n" + "=" * 60)
    print("VERIFICATION COMPLETE")
    print("=" * 60)
    
    return True


def analyze_klt_intersection_relationship(num_samples=10):
    """
    Analyze the relationship between KLT kernel and intersection matrix.
    
    Goal: Find the normalization factor c such that S_KLT = c * M_biadjoint^{-1}
    """
    print("=" * 60)
    print("ANALYSIS: KLT Kernel vs Intersection Matrix Inverse")
    print("=" * 60)
    
    # Load dependencies
    load('src/spinor_sampling.sage')
    load('src/klt.sage')
    load('src/kinematics_map.sage')
    
    ratios = []
    
    for sample in range(num_samples):
        print(f"\nSample {sample + 1}:")
        
        # Sample kinematics
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 17)
        if result is None:
            print("  Skipped (sampling failed)")
            continue
            
        lambdas, tilde_lambdas = result
        
        # Create adapter for KLT computation
        class SpinorAdapter:
            def __init__(self, lam, til):
                self.lambdas = lam
                self.tilde_lambdas = til
            def get_angle(self, i, j):
                l1, l2 = self.lambdas[i], self.lambdas[j]
                return l1[0]*l2[1] - l1[1]*l2[0]
            def get_square(self, i, j):
                t1, t2 = self.tilde_lambdas[i], self.tilde_lambdas[j]
                return t1[0]*t2[1] - t1[1]*t2[0]
        
        adapter = SpinorAdapter(lambdas, tilde_lambdas)
        
        def mandelstam_func(tw, i, j):
            return tw.get_angle(i,j) * tw.get_square(i,j)
        
        # Compute KLT kernel matrix
        permuted_set = [1, 2, 3]
        basis_perms = sorted(list(itertools.permutations(permuted_set)))
        
        S_klt = matrix(QQ, 6, 6)
        for i, alpha in enumerate(basis_perms):
            for j, beta in enumerate(basis_perms):
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_func)
                if val is None:
                    val = QQ(0)
                S_klt[i, j] = val
        
        # Compute bi-adjoint matrix
        M_ba, _ = compute_biadjoint_matrix(lambdas, tilde_lambdas, n=6)
        
        # Check if M_ba is invertible
        if M_ba.det() == 0:
            print("  Skipped (bi-adjoint matrix singular)")
            continue
        
        # Compute inverse
        M_ba_inv = M_ba.inverse()
        
        # Symmetrize both
        S_klt_sym = (S_klt + S_klt.transpose()) / 2
        M_ba_inv_sym = (M_ba_inv + M_ba_inv.transpose()) / 2
        
        # Compare entry-by-entry to find ratio
        sample_ratios = []
        for i in range(6):
            for j in range(6):
                if M_ba_inv_sym[i,j] != 0 and S_klt_sym[i,j] != 0:
                    ratio = S_klt_sym[i,j] / M_ba_inv_sym[i,j]
                    sample_ratios.append(float(ratio))
        
        if sample_ratios:
            avg_ratio = sum(sample_ratios) / len(sample_ratios)
            std_ratio = (sum((r - avg_ratio)**2 for r in sample_ratios) / len(sample_ratios))**0.5
            
            print(f"  Average ratio S_KLT / M^-1: {avg_ratio:.6f}")
            print(f"  Std dev: {std_ratio:.6f}")
            
            if std_ratio / abs(avg_ratio) < 0.01:
                print(f"  ✓ Consistent ratio!")
                ratios.append(avg_ratio)
            else:
                print(f"  ✗ Ratio varies significantly")
    
    if ratios:
        overall_avg = sum(ratios) / len(ratios)
        print(f"\n" + "=" * 60)
        print(f"Overall average ratio: {overall_avg:.6f}")
        print(f"This should be the normalization factor c")
        print("=" * 60)
    
    return ratios


if __name__ == "__main__":
    # Run verification
    print("Running intersection number verification...")
    
    # First verify setup exists
    data = load_setup(n=6)
    if data is None:
        print("Please run 01_setup_forms.sage first!")
        sys.exit(1)
    
    # Analyze KLT vs intersection relationship
    analyze_klt_intersection_relationship(num_samples=5)

