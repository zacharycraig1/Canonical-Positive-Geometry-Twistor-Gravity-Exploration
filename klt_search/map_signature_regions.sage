
from sage.all import *
from itertools import permutations
import numpy as np

# Load dependencies
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')

# Re-use KLT kernel for speed (as it matches M^-1 signature class)
# Or use M_biadjoint if KLT is still asymmetric/problematic?
# Previous results showed KLT symmetric part has (3,3).
# And M_biadjoint has (3,3).
# Let's use KLT symmetric part as it's faster than diagram summation?
# No, KLT kernel calculation is fast (sum over 3 terms * 3 terms?). 
# Actually KLT is faster.

def get_region_signature(adapter, mandelstam_func):
    # 1. Compute Signatures of Mandelstam Invariants (The "Region")
    # Planar set for n=6:
    # 2-particle: s12, s23, s34, s45, s56, s61
    # 3-particle: s123, s234, s345
    # (Note: s123 = s456, etc.)
    
    invariants = []
    names = []
    
    # s_i_i+1
    for i in range(6):
        j = (i + 1) % 6
        val = mandelstam_func(adapter, i, j)
        invariants.append(val)
        names.append(f"s{i+1}{j+1}")
        
    # s_i_i+1_i+2
    for i in range(3): # Only 3 independent ones? 
        # s123, s234, s345. 
        # s456 = s123.
        j = (i + 1) % 6
        k = (i + 2) % 6
        # s_ijk = s_ij + s_jk + s_ik? No.
        # s_ijk = (pi+pj+pk)^2 = s_ij + s_jk + s_ki
        s_ij = mandelstam_func(adapter, i, j)
        s_jk = mandelstam_func(adapter, j, k)
        s_ki = mandelstam_func(adapter, k, i) # s_13, etc.
        val = s_ij + s_jk + s_ki
        invariants.append(val)
        names.append(f"s{i+1}{j+1}{k+1}")
        
    # Region tuple: (+, -, +, ...)
    region_signs = tuple([1 if x > 0 else -1 for x in invariants])
    
    # 2. Compute KLT Kernel Signature
    # We use the standard KLT kernel function
    permuted_set = [1, 2, 3] # indices for legs 2,3,4
    basis_perms = sorted(list(permutations(permuted_set)))
    
    matrix_data = []
    for alpha in basis_perms:
        row = []
        for beta in basis_perms:
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_func)
            row.append(val)
        matrix_data.append(row)
        
    S = matrix(QQ, matrix_data)
    S_sym = (S + S.transpose()) / 2
    
    # Compute signature
    try:
        # Use RDF for eigenvalues
        evals = S_sym.change_ring(RDF).eigenvalues()
        n_pos = sum(1 for e in evals if e > 1e-6)
        n_neg = sum(1 for e in evals if e < -1e-6)
        n_zero = sum(1 for e in evals if abs(e) <= 1e-6)
        signature = (n_pos, n_neg, n_zero)
    except:
        signature = ("error", 0, 0)
        
    return region_signs, signature, names

def map_regions(num_samples=200):
    print(f"Mapping Kinematic Regions to Metric Signatures ({num_samples} samples)...")
    
    region_map = {} # Region -> {Signature -> count}
    
    for k in range(num_samples):
        try:
            # Use seed to ensure variety
            # sample_spinor_helicity_conserving uses random, so it's fine
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6)
            adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
            
            # Use simple mandelstam adapter
            def mandelstam_simple(tw, i, j):
                return tw.get_angle(i, j) * tw.get_square(i, j)
                
            signs, sig, names = get_region_signature(adapter, mandelstam_simple)
            
            # Store more data about region: s_ij values
            # We already have signs.
            
            if signs not in region_map:
                region_map[signs] = {}
            
            if sig not in region_map[signs]:
                region_map[signs][sig] = 0
            region_map[signs][sig] += 1
            
        except Exception as e:
            # print(e)
            continue
            
        if (k+1) % 50 == 0:
            print(f"Processed {k+1} samples...")
            
    # Report
    print(f"\nAnalysis of {len(region_map)} Regions:")
    
    # Sort by most frequent region
    sorted_regions = sorted(region_map.items(), key=lambda x: sum(x[1].values()), reverse=True)
    
    for signs, sig_counts in sorted_regions: # Show all regions to find Euclidean
        count = sum(sig_counts.values())
        if count < 2: continue # Skip rare ones
        
        print(f"\nRegion (count={count}):")
        # Print signs nicely
        sign_str = ""
        positive_vars = []
        negative_vars = []
        
        for i, s in enumerate(signs):
            name = names[i]
            if s > 0: positive_vars.append(name)
            else: negative_vars.append(name)
            
        print(f"  Positive: {', '.join(positive_vars)}")
        print(f"  Negative: {', '.join(negative_vars)}")
        print(f"  Signatures: {dict(sig_counts)}")
        
        # Identification
        # Euclidean: All s_ij < 0 (Physical Euclidean) or All > 0 (Anti-Euclidean)
        # Note: s_ijk signs also matter?
        # In Euclidean region, all planar s variables are negative?
        # Actually, in Euclidean kinematics, s_ij = -(p_i+p_j)^2 < 0.
        # So check if all 2-particle are negative.
        
        s2_signs = signs[:6] # first 6 are s_ij
        if all(s < 0 for s in s2_signs):
            print("  -> TYPE: Euclidean (All 2-pt negative)")
        elif all(s > 0 for s in s2_signs):
            print("  -> TYPE: Anti-Euclidean (All 2-pt positive)")
        else:
            # Check for Multi-Regge or specific scattering?
            # 2->4 scattering: s12 > 0, others negative?
            pass
            
class SpinorHelicityAdapter:
    def __init__(self, lambdas, tilde_lambdas):
        self.n = len(lambdas)
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas
    def get_angle(self, i, j):
        return ang_bracket(self.lambdas, i, j)
    def get_square(self, i, j):
        return sq_bracket(self.tilde_lambdas, i, j)

if __name__ == "__main__":
    map_regions(100)

