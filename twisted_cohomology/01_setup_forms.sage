#!/usr/bin/env sage
"""
Twisted Cohomology Setup for Scattering Amplitudes

This module sets up the twisted cohomology framework for computing
intersection numbers that define the KLT kernel for gravity amplitudes.

Key structures:
- Twist potential W = sum_{i<j} s_{ij} log(z_i - z_j)
- Connection 1-form omega = dW
- Parke-Taylor basis for twisted forms
- Scattering equations as critical points of W
"""
from sage.all import *
import itertools

def setup_twisted_forms(n=6):
    """
    Set up twisted cohomology forms for n-point scattering.
    
    Gauge fixing: (z_1, z_{n-1}, z_n) = (0, 1, infinity)
    Dynamic variables: z_2, ..., z_{n-2}
    
    Returns a dictionary with all necessary data for intersection computations.
    """
    print(f"=" * 60)
    print(f"Setting up Twisted Cohomology Forms for n={n}")
    print(f"=" * 60)
    
    # 1. Define Variables
    # Fix z_1=0, z_{n-1}=1, z_n=infinity
    # Variables: z_2, ..., z_{n-2}
    
    z_vars = {}
    z_symbolic = {}  # For symbolic computations
    var_names = []
    
    for i in range(2, n-1):
        name = f"z_{i}"
        z_vars[i] = var(name)
        z_symbolic[i] = z_vars[i]
        var_names.append(name)
        
    # Fixed values
    z_vars[1] = QQ(0)
    z_vars[n-1] = QQ(1)
    z_symbolic[1] = QQ(0)
    z_symbolic[n-1] = QQ(1)
    # z_vars[n] is infinity, treated separately
    
    num_dynamic = n - 3
    print(f"\nGauge fixing: z_1=0, z_{n-1}=1, z_{n}=∞")
    print(f"Dynamic variables ({num_dynamic}): {var_names}")
    
    # 2. Define Mandelstam Variables s_{ij}
    s_vars = {}
    s_var_list = []
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            name = f"s_{i}_{j}"
            s_vars[(i,j)] = var(name)
            s_vars[(j,i)] = s_vars[(i,j)]
            s_var_list.append((i, j))
            
    print(f"Mandelstam variables: {len(s_var_list)} independent s_{{ij}}")
    
    # 3. Define the Twist Potential W and its derivative (omega)
    # W = sum_{i<j} s_{ij} log(z_i - z_j)
    # omega_k = dW/dz_k = sum_{j!=k} s_{kj} / (z_k - z_j)
    # This is the scattering equation for particle k when set to zero
    
    omega_components = {}
    for k in range(2, n-1):  # Dynamic variables only
        val = SR(0)
        z_k = z_vars[k]
        
        # Sum over j from 1 to n-1 (excluding n=infinity and k itself)
        for j in range(1, n):
            if j == k: 
                continue
            z_j = z_vars[j]
            val += s_vars[(k,j)] / (z_k - z_j)
            
        # Term j = n (infinity) vanishes as 1/infinity -> 0
        # (Or equivalently, momentum conservation implies sum_j s_kj = 0)
        
        omega_components[k] = val
        
    print(f"\nTwist Connection omega = dW:")
    for k in sorted(omega_components.keys()):
        print(f"  omega_{k} = dW/dz_{k} = {omega_components[k]}")
        
    # 4. Define Parke-Taylor Basis Forms
    # Standard basis: permutations of dynamic indices {2, ..., n-2}
    # Full ordering: (1, alpha, n-1, n) with gauge-fixed positions
    # PT(alpha) = 1 / prod_{consecutive pairs} (z_u - z_v)
    # Terms with z_n (infinity) are dropped
    
    dynamic_indices = list(range(2, n-1))  # [2, 3, 4] for n=6
    perms = list(itertools.permutations(dynamic_indices))
    perms = sorted(perms)  # Canonical ordering
    
    pt_forms = {}
    pt_orderings = {}
    
    for alpha in perms:
        # Full cyclic ordering: 1, alpha..., n-1, n
        ordering = [1] + list(alpha) + [n-1, n]
        pt_orderings[alpha] = ordering
        
        # Compute denominator (product of consecutive differences)
        # Skip any pair involving z_n = infinity
        denom = SR(1)
        for i in range(len(ordering)):
            u = ordering[i]
            v = ordering[(i+1) % len(ordering)]
            
            # Skip terms with z_n (infinity)
            if u == n or v == n:
                continue
                
            val_u = z_vars[u]
            val_v = z_vars[v]
            denom *= (val_u - val_v)
            
        pt_forms[alpha] = 1 / denom
        
    print(f"\nParke-Taylor Basis ({len(perms)} forms):")
    for i, alpha in enumerate(perms):
        ordering_str = " -> ".join(map(str, pt_orderings[alpha]))
        print(f"  PT_{i} = PT{alpha}: ordering ({ordering_str})")
        if i < 2:  # Show first 2 explicit forms
            print(f"        = {pt_forms[alpha]}")
    
    # 5. Compute Scattering Equations (critical points of W)
    # These are omega_k = 0 for all k
    scattering_eqs = []
    for k in range(2, n-1):
        # Clear denominators to get polynomial equation
        eq = omega_components[k]
        # Multiply through by all denominators
        denom_prod = SR(1)
        for j in range(1, n):
            if j == k: 
                continue
            if j == n:
                continue
            denom_prod *= (z_vars[k] - z_vars[j])
        
        eq_poly = (eq * denom_prod).expand()
        scattering_eqs.append(eq_poly)
        
    print(f"\nScattering Equations (ω_k = 0):")
    for k, eq in zip(range(2, n-1), scattering_eqs):
        # Show simplified form
        print(f"  f_{k}(z) = 0")
    
    # 6. Compile all data
    data = {
        'n': n,
        'num_dynamic': num_dynamic,
        'dynamic_indices': dynamic_indices,
        'z_vars': z_vars,
        'z_symbolic': z_symbolic,
        'var_names': var_names,
        's_vars': s_vars,
        's_var_list': s_var_list,
        'omega': omega_components,
        'pt_forms': pt_forms,
        'pt_orderings': pt_orderings,
        'perms': perms,
        'scattering_eqs': scattering_eqs,
    }
    
    # Save to file
    save_path = f"twisted_cohomology/setup_n{n}.sobj"
    save(data, save_path)
    print(f"\nSaved to {save_path}")
    
    return data


def get_scattering_equation_solutions(data, s_values, precision=100):
    """
    Solve scattering equations numerically for given Mandelstam values.
    
    Args:
        data: Output from setup_twisted_forms
        s_values: Dictionary mapping (i,j) -> numerical value of s_ij
        precision: Number of bits for numerical computation
        
    Returns:
        List of solutions, each a dict mapping z_k -> value
    """
    n = data['n']
    z_vars = data['z_vars']
    omega = data['omega']
    
    # Substitute s values
    eqs = []
    for k in range(2, n-1):
        eq = omega[k]
        for (i,j), val in s_values.items():
            eq = eq.subs(data['s_vars'][(i,j)] == val)
        # Clear denominators
        denom_prod = 1
        for j in range(1, n):
            if j == k or j == n:
                continue
            denom_prod *= (z_vars[k] - z_vars[j])
        eq_poly = (eq * denom_prod).expand()
        eqs.append(eq_poly == 0)
    
    # Solve numerically
    z_var_list = [z_vars[k] for k in range(2, n-1)]
    
    try:
        # Use complex field for numerical solutions
        CC = ComplexField(precision)
        solutions = solve(eqs, z_var_list, solution_dict=True)
        return solutions
    except Exception as e:
        print(f"Solve failed: {e}")
        return []


def verify_setup_n6():
    """Verify the n=6 setup is correct."""
    print("\n" + "=" * 60)
    print("VERIFICATION: n=6 Setup")
    print("=" * 60)
    
    data = setup_twisted_forms(n=6)
    
    # Check dimensions
    assert data['n'] == 6, "n should be 6"
    assert data['num_dynamic'] == 3, "Should have 3 dynamic variables"
    assert len(data['perms']) == 6, "Should have 6 Parke-Taylor forms (3! = 6)"
    assert len(data['s_var_list']) == 15, "Should have 15 Mandelstam variables (6 choose 2)"
    
    print("\n✓ All dimension checks passed")
    
    # Verify Parke-Taylor forms are distinct
    pt_values = list(data['pt_forms'].values())
    assert len(set(str(pt) for pt in pt_values)) == 6, "PT forms should be distinct"
    
    print("✓ All Parke-Taylor forms are distinct")
    
    # Verify gauge fixing
    assert data['z_vars'][1] == 0, "z_1 should be 0"
    assert data['z_vars'][5] == 1, "z_5 should be 1"
    
    print("✓ Gauge fixing correct: z_1=0, z_5=1, z_6=∞")
    
    print("\n" + "=" * 60)
    print("n=6 SETUP COMPLETE")
    print("=" * 60)
    
    return data


if __name__ == "__main__":
    # Run both n=5 and n=6 setups
    print("Setting up n=5...")
    setup_twisted_forms(n=5)
    
    print("\n")
    
    print("Setting up n=6...")
    data = verify_setup_n6()









