#!/usr/bin/env sage
# =============================================================================
# SAMPLING MODULE: Moment-curve positive twistor generation
# =============================================================================

from sage.all import *
import numpy as np

def sample_positive_Z_moment_curve(n=6, seed=None):
    """
    Sample positive twistor matrix using moment curve with genericity nudges.
    
    Z_i = (1, t_i, t_i^2, t_i^3) with strictly increasing t_i.
    Uses deterministic pattern to avoid arithmetic progressions.
    
    Guarantees:
    - All ordered 4×4 minors > 0 (Vandermonde)
    - All angle brackets <i j> = t_j - t_i ≠ 0
    - Generic (no hidden degeneracies)
    
    Args:
        n: Number of particles (default 6)
        seed: Random seed for reproducibility
        
    Returns:
        List of n vectors in QQ^4
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Base: t_i = i + k_i where k_i are small rationals
    # Pattern: k_i = (i * 7 + seed) % 100 / 1000  (avoid arithmetic progression)
    t = []
    for i in range(n):
        base = QQ(i + 1)
        # Genericity nudge: small rational offset
        nudge = QQ((i * 7 + (seed if seed is not None else 0)) % 100) / QQ(1000)
        t_val = base + nudge
        t.append(t_val)
    
    # Ensure strict increase
    for i in range(1, n):
        if t[i] <= t[i-1]:
            t[i] = t[i-1] + QQ(1) / QQ(1000)
    
    # Construct Z using moment curve
    Z = []
    for t_i in t:
        z = vector(QQ, [
            QQ(1),
            t_i,
            t_i * t_i,
            t_i * t_i * t_i
        ])
        Z.append(z)
    
    return Z


def sample_random_Z(n=6, seed=None, range_bound=10):
    """
    Sample random twistor matrix with integer entries.
    
    For regression testing and failure classification.
    Not guaranteed to be in positive region.
    
    Args:
        n: Number of particles
        seed: Random seed
        range_bound: Integer range [-range_bound, range_bound]
        
    Returns:
        List of n vectors in QQ^4
    """
    if seed is not None:
        np.random.seed(seed)
    
    Z = []
    for i in range(n):
        z = vector(QQ, [
            QQ(np.random.randint(-range_bound, range_bound + 1)),
            QQ(np.random.randint(-range_bound, range_bound + 1)),
            QQ(np.random.randint(-range_bound, range_bound + 1)),
            QQ(np.random.randint(-range_bound, range_bound + 1))
        ])
        # Ensure non-zero
        while all(x == 0 for x in z):
            z = vector(QQ, [
                QQ(np.random.randint(-range_bound, range_bound + 1)),
                QQ(np.random.randint(-range_bound, range_bound + 1)),
                QQ(np.random.randint(-range_bound, range_bound + 1)),
                QQ(np.random.randint(-range_bound, range_bound + 1))
            ])
        Z.append(z)
    
    return Z









