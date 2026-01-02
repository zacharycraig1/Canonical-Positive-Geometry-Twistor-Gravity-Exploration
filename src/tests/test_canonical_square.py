import pytest
from sage.all import *
import sys
import os

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.canonical_form import CanonicalFormEvaluator

def test_square_triangulation_independence():
    """
    Test that a square evaluated with two different triangulations gives the same result.
    """
    # Square in P2 (z=1 plane)
    # V1=(1,1), V2=(-1,1), V3=(-1,-1), V4=(1,-1)
    # Homogeneous:
    Z1 = vector([1, 1, 1])
    Z2 = vector([-1, 1, 1])
    Z3 = vector([-1, -1, 1])
    Z4 = vector([1, -1, 1])
    
    # Triangulation 1: (Z1, Z2, Z3) + (Z1, Z3, Z4)
    # Note: Need consistent orientation (counter-clockwise)
    tris1 = [
        [Z1, Z2, Z3],
        [Z1, Z3, Z4]
    ]
    
    # Triangulation 2: (Z1, Z2, Z4) + (Z2, Z3, Z4)
    tris2 = [
        [Z1, Z2, Z4],
        [Z2, Z3, Z4]
    ]
    
    # Random W
    # W = vector([2, 3, 5]) 
    # This W hits a pole at (-1, -1, 1): -2 -3 + 5 = 0.
    # Use a safer W.
    W = vector([1, 2, 10]) # 1*x + 2*y + 10.
    # Checks:
    # (1,1,1) -> 1+2+10=13
    # (-1,1,1) -> -1+2+10=11
    # (-1,-1,1) -> -1-2+10=7
    # (1,-1,1) -> 1-2+10=9
    # All non-zero.
    
    val1 = 0
    for tri in tris1:
        val1 += CanonicalFormEvaluator.eval_simplex(tri, W)
        
    val2 = 0
    for tri in tris2:
        val2 += CanonicalFormEvaluator.eval_simplex(tri, W)
        
    # They should match
    assert abs(val1 - val2) < 1e-10
    
def test_polytope_eval_consistency():
    """
    Test that eval_polytope method matches manual triangulation sum.
    """
    # Square again
    points = [[1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1]]
    P = Polyhedron(vertices=points)
    
    W = vector([2, 3, 5])
    
    # Automatic evaluation
    # Note: eval_polytope sums terms. If Sage's triangulation is arbitrary, signs might mix.
    # We might need to take absolute values if the form is positive.
    # The canonical form of a convex polytope is positive for W inside the dual polytope.
    # W=[2,3,5] -> 2x + 3y + 5 > 0 ?
    # Vertices: 2(1)+3(1)+5=10 > 0.
    # 2(-1)+3(1)+5 = 6 > 0.
    # 2(-1)+3(-1)+5 = 0 (POLE!).
    # Need to avoid pole.
    
    W = vector([0, 0, 1]) # W.Z = 1 for all Z. Safe.
    
    val_auto = CanonicalFormEvaluator.eval_polytope(P, W)
    
    # Manual (Triangulation 1 is CCW)
    Z1 = vector([1, 1, 1])
    Z2 = vector([-1, 1, 1])
    Z3 = vector([-1, -1, 1])
    Z4 = vector([1, -1, 1])
    tris1 = [[Z1, Z2, Z3], [Z1, Z3, Z4]]
    
    val_manual = 0
    for tri in tris1:
        val_manual += CanonicalFormEvaluator.eval_simplex(tri, W)
        
    # Check if they match (possibly up to sign if auto triangulation is weird)
    assert abs(abs(val_auto) - abs(val_manual)) < 1e-10


