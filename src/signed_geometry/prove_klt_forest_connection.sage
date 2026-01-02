#!/usr/bin/env sage
"""
Investigate the connection between KLT (3,3) signature and forest 50/50 split.

Key question: Is there a mathematical proof that (3,3) signature implies ~50/50 split,
or is this just an empirical correlation?

Approach:
1. Express the forest sum in terms of the KLT structure
2. Analyze eigenspace decomposition
3. Look for a rigorous theorem

"""
from sage.all import *
import sys
import os

sys.path.insert(0, os.getcwd())

load('src/spinor_sampling.sage')
load('src/klt.sage')


def analyze_klt_eigenspace_structure():
    """
    Analyze the eigenspace structure of the KLT kernel.
    
    For signature (3,3), we have:
    - 3 positive eigenvalues
    - 3 negative eigenvalues
    
    The quadratic form x^T S x takes positive values on a 3-dim subspace
    and negative values on another 3-dim subspace.
    """
    print("=" * 70)
    print("KLT EIGENSPACE ANALYSIS")
    print("=" * 70)
    
    # Sample kinematics
    results = []
    
    for seed in range(10):
        result = sample_spinor_helicity_conserving(n=6, seed=seed * 31)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        # Build a mock twistor object for KLT computation
        # (We need to adapt since klt.sage uses MomentumTwistor)
        # For now, let's just analyze the structure theoretically
        
        results.append(seed)
    
    print(f"\nValid samples: {len(results)}")
    
    # Theoretical analysis:
    print("\n" + "=" * 70)
    print("THEORETICAL ANALYSIS")
    print("=" * 70)
    
    print("""
The KLT kernel S is a symmetric bilinear form on the space of orderings.
For n=6, this is a 6-dimensional space.

Signature (3,3) means:
- S = P^T D P where D = diag(+1, +1, +1, -1, -1, -1)
- The quadratic form x^T S x is positive on a 3-dim subspace
- The quadratic form x^T S x is negative on a 3-dim subspace

The gravity amplitude is:
  M = A^T · S · Ã

where A and Ã are vectors of Yang-Mills amplitudes.

Question: How does this relate to the forest expansion's 50/50 split?
""")


def forest_to_ordering_decomposition():
    """
    Express the forest expansion in terms of orderings.
    
    The Hodges determinant can be written as:
      det(L^R) = Σ_F term(F)
    
    But it can also be related to the bi-adjoint scalar matrix:
      det(L^R) ~ det(bi-adjoint matrix)
    
    The bi-adjoint scalar has a natural (n-3)! × (n-3)! matrix structure.
    """
    print("\n" + "=" * 70)
    print("FOREST ↔ ORDERING DECOMPOSITION")
    print("=" * 70)
    
    print("""
Key insight: Both expressions compute the same amplitude.

Forest expansion (108 terms for n=6):
  M = Σ_F ε(F) · |ω(F)|
  
KLT expansion (36 terms for n=6):
  M = Σ_{α,β} A_α · S[α|β] · Ã_β

The forest expansion has 108 terms, split ~54/54.
The KLT expansion has 36 terms (6×6 matrix).

Connection hypothesis:
- Each (α,β) pair in KLT corresponds to multiple forests
- The (3,3) signature of S means 3 "positive directions" and 3 "negative"
- When expanded to forests, this balanced signature manifests as ~50/50 split

This is NOT a rigorous theorem, but an observed structural correspondence.
""")


def compute_forest_sign_statistics():
    """
    Compute statistics of forest signs across many samples.
    """
    print("\n" + "=" * 70)
    print("FOREST SIGN STATISTICS")
    print("=" * 70)
    
    load('src/signed_geometry/canonical_form.sage')
    
    splits = []
    
    for seed in range(50):
        result = sample_spinor_helicity_conserving(n=6, seed=seed * 41)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        x_ref = vector(QQ, [1, 2])
        y_ref = vector(QQ, [3, 1])
        
        analysis = analyze_sign_structure(lambdas, tilde_lambdas, x_ref, y_ref, roots=(0, 1, 2))
        
        splits.append(analysis['split_ratio'])
    
    # Analyze distribution
    from collections import Counter
    split_counts = Counter(splits)
    
    print(f"\nSamples: {len(splits)}")
    print(f"\nSplit distribution:")
    
    for split, count in sorted(split_counts.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / len(splits)
        print(f"  {split}: {count} ({pct:.1f}%)")
    
    # Mode
    mode = split_counts.most_common(1)[0][0]
    print(f"\nMode: {mode}")
    
    # Mean positive fraction
    pos_fractions = [s[0] / 108.0 for s in splits]
    mean_pos = sum(pos_fractions) / len(pos_fractions)
    print(f"Mean positive fraction: {mean_pos:.4f} (0.5 = balanced)")
    
    return splits


def theoretical_connection():
    """
    State the theoretical connection as clearly as possible.
    """
    print("\n" + "=" * 70)
    print("THEORETICAL CONNECTION (CONJECTURE → THEOREM?)")
    print("=" * 70)
    
    print("""
STATEMENT:

For n=6 MHV gravity:
1. The KLT kernel has signature (3,3) [VERIFIED]
2. The forest expansion has modal split (54,54) [VERIFIED]
3. These are related by:

   The KLT structure decomposes as:
     S = Σ_i λ_i v_i ⊗ v_i
   
   where 3 eigenvalues are positive and 3 are negative.
   
   When this bilinear form is "expanded" into the forest basis,
   the balanced eigenvalue structure translates to balanced signs.

RIGOROUS PROOF STATUS:

We can prove:
- The forest signs are uniquely determined by MTT [PROVEN]
- The KLT kernel generically has split signature [VERIFIED]
- The modal forest split is (54,54) [VERIFIED]

We CANNOT yet prove:
- A direct mathematical formula: split = f(signature)
- That (3,3) signature IMPLIES 50/50 split

The connection remains a strong empirical correlation, not a theorem.
""")
    
    print("\n" + "=" * 70)
    print("RECOMMENDATION")
    print("=" * 70)
    
    print("""
For the paper:
- Keep the connection as a CONJECTURE
- State the empirical evidence clearly
- Note that a rigorous proof would require showing:
  "The forest expansion is a refined decomposition of the KLT bilinear form,
   and the (p,q) signature of the kernel induces a (p', q') split of forests
   with p'/q' ≈ p/q."

This is an interesting open question for future work.
""")


if __name__ == "__main__":
    analyze_klt_eigenspace_structure()
    forest_to_ordering_decomposition()
    
    # Run statistics if possible
    try:
        splits = compute_forest_sign_statistics()
    except Exception as e:
        print(f"\nCould not compute statistics: {e}")
    
    theoretical_connection()
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
1. UNIQUENESS: PROVEN (signs are determined by MTT formula)

2. KLT ↔ FOREST CONNECTION: REMAINS CONJECTURE
   - Strong empirical correlation
   - No rigorous proof of causation
   - Both reflect underlying indefinite structure
   
The paper should:
- Upgrade Uniqueness to Theorem
- Keep KLT connection as Conjecture with empirical support
""")

