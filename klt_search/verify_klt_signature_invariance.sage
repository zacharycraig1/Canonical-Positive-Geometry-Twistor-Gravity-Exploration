
from sage.all import *
from itertools import permutations
import os

# Load necessary sage files directly
# Paths are relative to the workspace root when run from root
load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')

class SpinorHelicityAdapter:
    """
    Adapts lambdas/tilde_lambdas to the interface expected by KLT.
    """
    def __init__(self, lambdas, tilde_lambdas):
        self.n = len(lambdas)
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas

    def get_angle(self, i, j):
        return ang_bracket(self.lambdas, i, j)

    def get_square(self, i, j):
        return sq_bracket(self.tilde_lambdas, i, j)

def mandelstam_adapter(twistor, i, j):
    return twistor.get_angle(i, j) * twistor.get_square(i, j)

def check_signature_invariance(num_samples=50):
    print(f"Checking signature invariance across {num_samples} samples...")
    
    signatures = {} # signature -> count
    regions = set() # Set of sign tuples
    
    permuted_set = [1, 2, 3] # These are indices in 0-based system corresponding to legs 2,3,4
    # The KLT kernel is indexed by permutations of these legs.
    # Wait, in klt.sage:
    #   alpha, beta are permutations of {1,2,3} (0-based indices)
    #   Fixed legs: {0, 4, 5}
    
    basis_perms = sorted(list(permutations(permuted_set)))
    n_basis = len(basis_perms)
    
    for k in range(num_samples):
        # Sample kinematics
        # Use a new seed each time implicitly or explicitly
        # sample_spinor_helicity_conserving uses random module
        try:
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6)
        except Exception:
            # Sometime sampling fails
            continue
            
        adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
        
        # 1. Compute Matrix
        matrix_data = []
        for i, alpha in enumerate(basis_perms):
            row = []
            for j, beta in enumerate(basis_perms):
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_adapter)
                row.append(val)
            matrix_data.append(row)
            
        S = matrix(QQ, matrix_data)
        
        # Symmetrize S
        S_sym = (S + S.transpose()) / 2
        
        # 2. Compute Signature of Symmetric Part
        # Eigenvalues of real symmetric matrix are real.
        
        try:
            evals = S_sym.eigenvalues()
        except Exception as e:
            # Fallback to RDF
            S_float = S_sym.change_ring(RDF)
            evals = S_float.eigenvalues()
            
        n_pos = 0
        n_neg = 0
        n_zero = 0
        n_complex = 0 # Should be 0 for symmetric matrix
        
        for e in evals:
            if e.imag() > 1e-10: # Check for complex
                n_complex += 1
            elif e.real() > 1e-10:
                n_pos += 1
            elif e.real() < -1e-10:
                n_neg += 1
            else:
                n_zero += 1
                
        sig = (n_pos, n_neg, n_zero, n_complex)
        signatures[sig] = signatures.get(sig, 0) + 1
        
        # 3. Identify Region
        # Planar Mandelstams: s12, s23, s34, s45, s56, s61, s123, s234, s345
        mandelstams = []
        # Two particle
        pairs = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
        for i, j in pairs:
            mandelstams.append(mandelstam_adapter(adapter, i, j))
        
        # Three particle
        triples = [(0,1,2), (1,2,3), (2,3,4)]
        for i, j, l in triples:
            # s_ijl = s_ij + s_il + s_jl
            s = mandelstam_adapter(adapter, i, j) + \
                mandelstam_adapter(adapter, i, l) + \
                mandelstam_adapter(adapter, j, l)
            mandelstams.append(s)
            
        signs = tuple([1 if x > 0 else -1 for x in mandelstams])
        regions.add(signs)
        
        if (k+1) % 10 == 0:
            print(f"Processed {k+1} samples. Regions found: {len(regions)}")

    print("\nResults:")
    print(f"Total samples: {num_samples}")
    print(f"Total unique kinematic regions visited: {len(regions)}")
    print("Signatures (pos, neg, zero, complex):")
    for sig, count in signatures.items():
        print(f"  {sig}: {count}")

    # Check basis independence briefly (as mentioned in next steps)
    # We can do this in another script, but if we see issues here...
    
if __name__ == "__main__":
    check_signature_invariance(100)

