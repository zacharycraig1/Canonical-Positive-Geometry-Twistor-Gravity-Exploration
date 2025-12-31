#!/usr/bin/env sage
from sage.all import *
load('src/klt.sage')
load('src/hodges.sage') # for MomentumTwistor
load('src/spinor_helicity.sage') # for physical mandelstam
load('src/kinematics_map.sage')

def analyze_kernel():
    print("Analyzing KLT Kernel Matrix Properties (n=6)")
    print("==========================================")
    
    # Setup
    # Permutations of {2,3,4} (indices 1,2,3). Fixed {1,5,6} (indices 0,4,5).
    import itertools
    perm_indices = [1, 2, 3]
    perms = sorted(list(itertools.permutations(perm_indices)))
    dim = len(perms)
    print(f"Kernel Dimension: {dim}x{dim}")
    
    # Generate Physical Point
    twistor = MomentumTwistor(n=6, check_domain=False)
    lambdas, tilde_lambdas, _ = extract_spinors_from_twistor(twistor)
    if lambdas is None: 
        print("Failed to extract spinors")
        return
        
    adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
    
    def mandelstam_func(tw, i, j):
        # Use adapter closure or pass adapter as tw
        # The klt_momentum_kernel_6pt expects 'twistor' as third arg
        # and 'mandelstam_func' as fourth. 
        # func(tw, i, j).
        # We will wrap it.
        return adapter.get_angle(i,j) * adapter.get_square(i,j)

    # Compute Matrix S
    S = matrix(QQ, dim, dim)
    
    for i in range(dim):
        for j in range(dim):
            alpha = list(perms[i])
            beta = list(perms[j])
            
            # Note: klt_momentum_kernel_6pt implementation in src/klt.sage
            # uses indices 0-based. 
            # Our perms are [1,2,3] (indices for particles 2,3,4).
            # The function expects perm indices.
            
            val = klt_momentum_kernel_6pt(alpha, beta, adapter, mandelstam_func)
            if val is None:
                print(f"Failed to compute S[{i},{j}]")
                return
            S[i,j] = val
            
    print("\nKernel Matrix S (First 3x3):")
    print(S[:3,:3])
    
    # Properties
    is_sym = S.is_symmetric()
    print(f"\nSymmetric: {is_sym}")
    
    if not is_sym:
        print("Diff:", (S - S.transpose()).norm())
        print("Analyzing Symmetric Part S_sym = 0.5 * (S + S.T)")
        S_sym = 0.5 * (S + S.transpose())
    else:
        S_sym = S
        
    # Rank
    print(f"Rank of S: {S.rank()}")
    print(f"Rank of S_sym: {S_sym.rank()}")
    
    # Eigenvalues (if symmetric, or just in general)
    try:
        # Convert to RDF for numerical eigenvalues
        S_float = S_sym.change_ring(RDF)
        eigs = S_float.eigenvalues()
        print("\nEigenvalues of S_sym:")
        for e in eigs:
            print(f"  {e:.4e}")
            
        # Check Definiteness
        # Check signs (ignore small numerical noise around 0)
        pos = sum(1 for e in eigs if e > 1e-6)
        neg = sum(1 for e in eigs if e < -1e-6)
        zeros = sum(1 for e in eigs if abs(e) <= 1e-6)
        
        print(f"\nSignature (pos, neg, zero): ({pos}, {neg}, {zeros})")
            
    except Exception as e:
        print(f"Eigenvalue error: {e}")

if __name__ == "__main__":
    analyze_kernel()

