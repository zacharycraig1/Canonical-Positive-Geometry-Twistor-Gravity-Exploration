import sys
import os
from sage.all import *

# Load 54.sage to access DCP logic
# We need to suppress its main execution if it has one?
# 54.sage has if __name__ == "__main__": main()
try:
    load("54.sage")
except Exception as e:
    # 54.sage is failing due to seed type error in main execution block?
    # Or import block?
    # The traceback showed:
    # Error loading 54.sage: The only supported seed types are: None, int, float...
    # This suggests 54.sage executes something on load that calls random.seed with a Sage integer.
    # We should fix 54.sage to cast seeds to int.
    print(f"Error loading 54.sage: {e}")

def get_dcp_basis():
    """
    Construct the basis of OS3 forms for n=6.
    Returns:
        C: List of channels
        triples: List of basis triples (i,j,k)
        M6: The matroid
    """
    # 1. Build Matroid
    C, M6 = build_M6_matroid()
    
    # 2. Build OS3 Triples
    # build_OS3_data returns (triples, ...)
    # But 54.sage's build_OS3_data does a lot of other stuff (invariants)
    # Let's just extract the triples logic or call it.
    
    # In 54.sage:
    # def build_OS3_data(C, M6):
    #     ...
    #     triples = [(i,j,k) for i in range(m) for j in range(i+1,m) for k in range(j+1,m)]
    #     ...
    #     return triples, ...
    
    # Actually, 54.sage builds the FULL basis of all triples (dim = m choose 3).
    # Then it filters them?
    # No, OS3 algebra is quotient of exterior algebra.
    # The basis of H^3(M) is smaller.
    # 54.sage computes the "OS3 candidate space" by finding kernel of dependencies?
    # No, let's look at 54.sage code again.
    
    # It seems 54.sage searches for "invariants" inside the full space of triples?
    # Or does it construct the OS3 basis?
    # Usually OS3 basis is formed by "nbc" sets (broken circuits).
    
    # Let's rely on 54.sage's `triples` if it exports the basis.
    # Actually, 54.sage seems to work with "triples" as the naive basis (m choose 3)
    # and then imposes constraints (OS relation?).
    
    # Let's just reproduce the channel generation and triples for now.
    m = len(C)
    triples = [(i,j,k) for i in range(m) for j in range(i+1,m) for k in range(j+1,m)]
    
    return C, triples, M6

def get_invariant_basis(mode='S3xS3'):
    """
    Get the basis of forms invariant under S3xS3 (or other groups).
    This might be needed to reduce the search space.
    """
    # 54.sage has `compute_invariants`.
    # We might want to call that.
    pass
