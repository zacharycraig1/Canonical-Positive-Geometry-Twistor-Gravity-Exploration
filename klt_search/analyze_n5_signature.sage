#!/usr/bin/env sage
from sage.all import *
load('src/klt.sage')

def analyze_kernel_n5():
    print("Analyzing KLT Kernel Matrix Properties (n=5)")
    print("==========================================")
    
    # n=5
    # Perms of {2,3} (indices 1,2). Fixed {1,4,5} (indices 0,3,4).
    # Wait, klt module is hardcoded for 6pt in some places?
    # No, klt_momentum_kernel_6pt is specific. 
    # We need a general klt kernel function or verify if existing one works.
    # src/klt.sage: klt_momentum_kernel_6pt is hardcoded loop range(3).
    
    # We need to implement generic KLT kernel for this check.
    
    def klt_kernel_generic(alpha, beta, twistor, mandelstam_func):
        # alpha, beta are permutations of {2, ..., n-2}
        # Fixed 1, n-1, n.
        # This implementation matches the n=6 one but generic
        
        # Determine permuted set from alpha
        perm_set = sorted(list(alpha))
        k = len(alpha) # n-3
        
        # Build position map for beta
        pos_in_beta = {val: idx for idx, val in enumerate(beta)}
        
        def theta_beta(a, b):
            return 1 if pos_in_beta[a] > pos_in_beta[b] else 0
            
        kernel = QQ(1)
        
        # Formula: prod_{i=1}^k (s_{1, alpha[i]} + sum_{j<i} theta(alpha[j], alpha[i]) s_{alpha[j], alpha[i]})
        # Indices in paper usually 1-based.
        # Our alpha is 0-based indices from the code.
        # Fixed leg 1 is index 0.
        
        fixed_leg_1 = 0
        
        for i in range(k):
            # Term corresponding to alpha[i]
            # s_{1, alpha[i]}
            term = mandelstam_func(twistor, fixed_leg_1, alpha[i])
            
            # Sum over j < i
            for j in range(i):
                if theta_beta(alpha[j], alpha[i]):
                    s_ij = mandelstam_func(twistor, alpha[j], alpha[i])
                    term += s_ij
            
            kernel *= term
            
        return kernel

    # Setup n=5
    import itertools
    perm_indices = [1, 2] # 2 particles permuted
    perms = sorted(list(itertools.permutations(perm_indices)))
    dim = len(perms)
    print(f"Dimension: {dim}x{dim}") # Should be 2x2
    
    # Mock Mandelstams (random symmetric matrix for s_ij)
    # s_ij variables.
    n = 5
    s_vals = {}
    for i in range(n):
        for j in range(i+1, n):
            import random
            val = QQ(random.randint(1, 100))
            s_vals[(i,j)] = val
            s_vals[(j,i)] = val
            
    def m_func(tw, i, j):
        return s_vals.get((i,j), 0)
        
    S = matrix(QQ, dim, dim)
    for i in range(dim):
        for j in range(dim):
            S[i,j] = klt_kernel_generic(perms[i], perms[j], None, m_func)
            
    print("Matrix S (n=5):")
    print(S)
    
    S_sym = 0.5 * (S + S.transpose())
    eigs = S_sym.eigenvalues()
    print("Eigenvalues:", eigs)
    
    # Check signature
    pos = sum(1 for e in eigs if e > 0)
    neg = sum(1 for e in eigs if e < 0)
    print(f"Signature: ({pos}, {neg})")

if __name__ == "__main__":
    analyze_kernel_n5()









