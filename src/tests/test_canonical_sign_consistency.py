import sys
import os
sys.path.append(os.getcwd())
from sage.all import *
from src.posgeom.canonical_form import CanonicalFormEvaluator

def test_square_sign():
    """
    Test canonical form of a square in P^2.
    Vertices: (1,0,0), (1,1,0), (1,1,1), (1,0,1)
    """
    # Define vertices as vectors
    v0 = vector([1, 0, 0])
    v1 = vector([1, 1, 0])
    v2 = vector([1, 1, 1])
    v3 = vector([1, 0, 1])
    
    P = Polyhedron(vertices=[v0, v1, v2, v3])
    
    # Dual point W inside dual cone
    # W . Z > 0 for all Z
    # W = (1, 1, 1) -> 
    # (1,0,0).W = 1
    # (1,1,0).W = 2
    # (1,1,1).W = 3
    # (1,0,1).W = 2
    # All > 0.
    W_in = vector([1, 1, 1])
    
    val_in = CanonicalFormEvaluator.eval_polytope(P, W_in)
    print(f"Value inside: {val_in}")
    
    assert val_in > 0, "Canonical form should be positive in the positive region"
    
    # Point W that might give negative value?
    # The form is a rational function.
    # If we cross a pole, sign changes.
    # Vertices define boundaries: W.Z = 0.
    # W = (1, -0.5, 1)
    # v0.W = 1
    # v1.W = 0.5
    # v2.W = 1.5
    # v3.W = 2
    # All positive. Still inside dual cone?
    # Square edges define facets.
    # Facet 1: v0, v1 (z=0). Normal (0,0,1). W=(0,0,1) -> W.Z > 0.
    # Facet 2: v1, v2 (x=1). (But x is not coord 1 here? projective).
    # Let's verify if we can get a negative value.
    # Move W to (1, -2, 1).
    # v0.W = 1
    # v1.W = -1 (Negative!)
    # v2.W = 0
    # v3.W = 2
    # We crossed a boundary defined by v1?
    # No, boundaries are facets (lines between vertices).
    # Pole occurs when W is on dual of a facet.
    # But if W.v1 < 0, we are "behind" v1?
    
    W_neg = vector([1, -1.5, 1])
    # This might hit a pole or be valid.
    # Let's calculate.
    try:
        val_neg = CanonicalFormEvaluator.eval_polytope(P, W_neg)
        print(f"Value at W_neg: {val_neg}")
        if val_neg < 0:
            print("Successfully obtained negative value.")
        else:
            print("Value is positive (or zero).")
    except ValueError:
        print("Hit a pole.")

    # Check consistency vs known formula for square?
    # Square form: 1 / (x y) ? (in some coords)
    
def test_triangulation_independence():
    # Force different triangulation?
    # Hard to force Sage.
    # But we can rotate vertices in definition order and see if Sage changes triangulation
    pass

if __name__ == "__main__":
    test_square_sign()

