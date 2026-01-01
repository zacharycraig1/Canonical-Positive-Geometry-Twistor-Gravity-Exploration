import pytest
from sage.all import *
import sys
import os

# Ensure src is in path
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.posgeom.canonical_form import CanonicalFormEvaluator

def test_simplex_formula():
    """
    Test that eval_simplex matches the closed form formula.
    """
    # Define a standard 2-simplex in P2
    Z0 = vector([1, 0, 0])
    Z1 = vector([0, 1, 0])
    Z2 = vector([0, 0, 1])
    vertices = [Z0, Z1, Z2]
    
    # Dual point W
    W = vector([1, 1, 1])
    
    # Expected: det(I) / (1*1*1) = 1
    val = CanonicalFormEvaluator.eval_simplex(vertices, W)
    assert val == 1
    
    # Scale W
    W2 = vector([2, 2, 2])
    # Expected: 1 / (2*2*2) = 1/8
    # BUT wait, the form is homogenous of degree -d-1?
    # Dim = 2. d+1 = 3. Yes.
    val2 = CanonicalFormEvaluator.eval_simplex(vertices, W2)
    assert val2 == 1/8

def test_simplex_orientation():
    """
    Test that swapping vertices flips the sign of eval_simplex.
    """
    Z0 = vector([1, 0, 0])
    Z1 = vector([0, 1, 0])
    Z2 = vector([0, 0, 1])
    vertices = [Z0, Z1, Z2]
    vertices_swapped = [Z1, Z0, Z2]
    
    W = vector([1, 1, 1])
    
    val1 = CanonicalFormEvaluator.eval_simplex(vertices, W)
    val2 = CanonicalFormEvaluator.eval_simplex(vertices_swapped, W)
    
    assert val1 == -val2



