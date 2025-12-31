#!/usr/bin/env sage
from sage.all import *

def setup_twisted_forms(n=5):
    print(f"Setting up Twisted Cohomology Forms for n={n}")
    
    # 1. Define Variables
    # Fix z_1=0, z_{n-1}=1, z_n=infinity
    # Variables: z_2, ..., z_{n-2}
    # Indices in Python: z[2]...z[n-2]
    
    # We'll use a dictionary for z variables
    z_vars = {}
    var_names = []
    for i in range(2, n-1):
        name = f"z_{i}"
        z_vars[i] = var(name)
        var_names.append(name)
        
    # Fixed values
    z_vars[1] = 0
    z_vars[n-1] = 1
    # z_vars[n] is infinity, treated separately in forms
    
    print(f"Variables: {var_names}")
    
    # 2. Define Mandelstam Variables s_{ij}
    # For Twisted Cohomology, we need symmetric s_{ij}
    s_vars = {}
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            name = f"s_{i}_{j}"
            s_vars[(i,j)] = var(name)
            s_vars[(j,i)] = s_vars[(i,j)]
            
    # 3. Define the Twist \omega
    # \omega = d \log(\Phi) = \sum_{1<=i<j<=n} s_{ij} d \log(z_i - z_j)
    # \omega = \sum_{k} \omega_k dz_k
    # We only sum over dynamic variables z_k
    
    # Note: terms involving z_n (infinity) behave as -d \log(z_i) ? 
    # Or using Mobius invariance condition sum s_{ij} = 0.
    # Potential function W = sum_{i<j} s_{ij} log(z_i - z_j)
    # dW = sum_k (dW/dz_k) dz_k
    
    # Derivatives w.r.t dynamic z_k:
    # dW/dz_k = sum_{j!=k} s_{kj} / (z_k - z_j)
    
    omega_components = {}
    for k in range(2, n-1): # Dynamic variables
        val = 0
        z_k = z_vars[k]
        
        # Terms j != n
        for j in range(1, n):
            if j == k: continue
            z_j = z_vars[j]
            val += s_vars[(k,j)] / (z_k - z_j)
            
        # Term j = n (infinity)
        # Using momentum conservation sum_{j} s_{kj} = 0
        # The term is implicitly handled or vanishes? 
        # Actually, standard formula omits z_n terms if sum s_{kj} = 0
        # But let's keep it explicit if needed. Usually limit z_n -> inf implies term vanishes (1/inf).
        # So we just sum j from 1 to n-1.
        
        omega_components[k] = val
        
    print("\nTwist Components (omega):")
    for k in omega_components:
        print(f"  dz_{k}: {omega_components[k]}")
        
    # 4. Define Parke-Taylor Forms
    # PT(alpha) = 1 / (z_{a1}-z_{a2})(z_{a2}-z_{a3})...(z_{an}-z_{a1})
    # We need to account for fixed points.
    # PT form is usually defined on the moduli space.
    # Standard basis: alpha is perm of {2, ..., n-2} (dynamic points)
    # Ordering (1, alpha, n-1, n)
    
    import itertools
    dynamic_indices = list(range(2, n-1))
    perms = list(itertools.permutations(dynamic_indices))
    
    pt_forms = {}
    
    for alpha in perms:
        # Full ordering: 1, alpha..., n-1, n
        ordering = [1] + list(alpha) + [n-1, n]
        
        # Denominator product
        denom = 1
        for i in range(len(ordering)):
            u = ordering[i]
            v = ordering[(i+1)%len(ordering)]
            
            # If one is n (infinity), the term (z_u - z_v) becomes 1 (or is removed from PT factor)
            # Standard PT factor rule: drop terms with infinity?
            # Or Faddeev-Popov determinant factor?
            # Standard rule: PT(1,2,...,n) = 1/((z1-z2)...(zn-z1))
            # With fixed z1, z(n-1), zn:
            # The measure is d^{n-3}z.
            # The PT factor reduces to 1/((z1-z2)...(z(n-1)-1)) roughly.
            # Let's use the explicit product but treat (z - z_n) as 1?
            
            # Actually, proper way: 
            # PT_factor = 1 / prod_{i} (z_i - z_{i+1})
            # Drop terms with z_n.
            # And multiply by (z_1 - z_{n-1})(z_{n-1} - z_n)(z_n - z_1) -> (z_1 - z_{n-1})
            # To compensate for gauge fixing?
            
            if u == n or v == n:
                continue
                
            val_u = z_vars[u]
            val_v = z_vars[v]
            denom *= (val_u - val_v)
            
        pt_forms[alpha] = 1/denom
        
    print(f"\nParke-Taylor Forms ({len(perms)} basis elements):")
    for alpha in perms[:3]: # Show first 3
        print(f"  PT{alpha}: {pt_forms[alpha]}")

    # Save for next steps
    data = {
        'n': n,
        'omega': omega_components,
        'pt_forms': pt_forms,
        'z_vars': var_names,
        's_vars': [(k[0], k[1]) for k in s_vars.keys()]
    }
    save(data, f"twisted_cohomology/setup_n{n}.sobj")

if __name__ == "__main__":
    setup_twisted_forms(n=5)





