#!/usr/bin/env sage
from sage.all import *
load('src/klt.sage')

def analyze_klt_kernel_structure():
    print("Analyzing KLT Kernel Structure for n=6")
    print("=======================================")
    
    # 1. Define Symbolic Variables
    # We need s_ij for the kernel.
    # Kernel uses s_{0,alpha[i]} and s_{alpha[j],alpha[i]}
    # Indices are 0,1,2,3,4,5.
    
    s_vars = {}
    for i in range(6):
        for j in range(i+1, 6):
            s_vars[(i,j)] = var(f"s_{i+1}{j+1}")
            s_vars[(j,i)] = s_vars[(i,j)]
            
    def symbolic_mandelstam(twistor, i, j):
        return s_vars.get((i,j), 0)
        
    # 2. Iterate over Basis Pairs
    # Permuted set {1,2,3} (indices 1,2,3). Fixed 0,4,5.
    import itertools
    perms = sorted(list(itertools.permutations([1, 2, 3])))
    
    print(f"Basis size: {len(perms)}")
    
    # We will analyze diagonal and off-diagonal kernel elements
    
    for i_a, alpha in enumerate(perms):
        for i_b, beta in enumerate(perms):
            # Compute Kernel S[alpha|beta]
            S_sym = klt_momentum_kernel_6pt(list(alpha), list(beta), None, symbolic_mandelstam)
            
            # Simplify/Expand
            S_expanded = S_sym.expand()
            
            if S_expanded != 0:
                print(f"\nKernel S[{alpha} | {beta}]:")
                print(f"  {S_expanded}")
                
                # Check degree
                # It should be degree 3 in s
                # We can check by inspection of the string or polynomial degree
                pass

if __name__ == "__main__":
    analyze_klt_kernel_structure()

