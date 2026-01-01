
import sys
import os
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.canonical_form import CanonicalFormEvaluator
from src.posgeom.forest_polytope import get_forest_exponents

def test_triangulation_invariance():
    print("S4: Testing Triangulation Invariance & Scaling...")
    
    # Setup Forest Polytope (n=5, roots={0,1})
    n = 5
    roots = [0, 1]
    exponents, _ = get_forest_exponents(n, roots)
    
    # Vertices in R^10 (since 5 choose 2 = 10 edges)
    # Dimension should be 5 - 2 - 1 = 2 ?? No.
    # Vertices of forest polytope are indicator vectors.
    P = Polyhedron(vertices=exponents)
    print(f"Polytope dim: {P.dim()}, Ambient dim: {len(exponents[0])}")
    
    # Generate random W in dual space (Ambient dim + 1)
    dim_ambient = len(exponents[0])
    
    # Run 5 trials
    for i in range(5):
        set_random_seed(i + 100) # New seed
        
        # Try to find a valid W
        W = None
        Z_vecs = [vector(QQ, [1] + list(v)) for v in exponents]
        
        for attempt in range(100):
            # Use larger range and maybe rational denominator to avoid accidental zeros
            cand = vector(QQ, [QQ(randint(-100, 100))/randint(1,5) for _ in range(dim_ambient + 1)])
            if all(cand.dot_product(z) != 0 for z in Z_vecs):
                W = cand
                break
        
        if W is None:
            print(f"Could not find non-singular W for trial {i}")
            continue
            
        print(f"\nTrial {i}: W = {W[:3]}...")
        
        # 1. Base Evaluation
        try:
            val1 = CanonicalFormEvaluator.eval_polytope(P, W)
            print(f"  Val (Order 1): {val1}")
        except Exception as e:
            print(f"  Error Order 1: {e}")
            continue
            
        # 2. Permuted Vertex Order (should induce different triangulation)
        import random
        shuffled_exponents = list(exponents)
        random.shuffle(shuffled_exponents)
        P2 = Polyhedron(vertices=shuffled_exponents)
        
        try:
            val2 = CanonicalFormEvaluator.eval_polytope(P2, W)
            print(f"  Val (Order 2): {val2}")
        except Exception as e:
            print(f"  Error Order 2: {e}")
            val2 = None

        if val1 is not None and val2 is not None:
            diff = val1 - val2
            if abs(diff) < 1e-10:
                print("  [PASS] Triangulation Invariance")
            else:
                print(f"  [FAIL] Difference: {diff}")
                
        # 3. Scaling Test
        # Omega(cW) = c^{-(d+1)} Omega(W) where d = P.dim()
        # Note: eval_polytope projects P to R^d.
        # So effective homogeneous degree is -(d+1).
        
        scale = QQ(2)
        W_scaled = scale * W
        
        val_scaled = CanonicalFormEvaluator.eval_polytope(P, W_scaled)
        
        d = P.dim()
        expected_ratio = scale ** (-(d + 1))
        
        ratio = val_scaled / val1
        print(f"  Scaling Ratio: {ratio} (Expected 2^{-d-1} = {expected_ratio})")
        
        if abs(ratio - expected_ratio) < 1e-10:
            print("  [PASS] Scaling Property")
        else:
            print("  [FAIL] Scaling Property")

if __name__ == "__main__":
    test_triangulation_invariance()
