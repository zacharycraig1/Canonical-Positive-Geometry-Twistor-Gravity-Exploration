#!/usr/bin/env sage
"""
BCFW Cell Enumeration and Canonical Forms
==========================================

This module implements BCFW cells for MHV gravity amplitudes.

For MHV gravity (k=2), BCFW cells correspond to triangulations of the n-gon.
Each cell contributes a d-log form, and the sum should equal Hodges.

Key insight from KLT:
    M_gravity = ∑_{perms} S[α|β] A_YM(α) A_YM(β)

This suggests gravity lives in a PRODUCT space, related to Gr(2,n) × Gr(2,n).
"""

from sage.all import *
from itertools import combinations, permutations
import sys
import os

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from src.amplituhedron.momentum_twistor import MomentumTwistorData
except ImportError:
    # Handle load order issues in Sage
    pass


def catalan_number(n):
    """Compute n-th Catalan number."""
    return binomial(2*n, n) // (n + 1)


def triangulations_of_polygon(n):
    """
    Enumerate all triangulations of a convex n-gon.
    
    A triangulation divides the n-gon into (n-2) triangles using (n-3) diagonals.
    The number of triangulations is Catalan(n-2).
    
    Returns list of triangulations, each as a list of diagonals (i, j) with i < j.
    """
    if n < 3:
        return []
    if n == 3:
        return [[]]  # Triangle needs no diagonals
    if n == 4:
        # Two triangulations: (0,2) or (1,3)
        return [[(0, 2)], [(1, 3)]]
    
    # Recursive construction for n ≥ 5
    # Fix edge (0, n-1), then vertex v is connected to both 0 and n-1
    triangulations = []
    
    for v in range(1, n - 1):
        # Triangulate left polygon: 0, 1, ..., v
        left = triangulations_of_polygon(v + 1) if v + 1 >= 3 else [[]]
        
        # Triangulate right polygon: v, v+1, ..., n-1
        # Relabel vertices: v→0, v+1→1, ..., n-1→(n-1-v)
        right_size = n - v
        right_raw = triangulations_of_polygon(right_size) if right_size >= 3 else [[]]
        
        # Relabel right diagonals back to original indices
        right = []
        for tri in right_raw:
            relabeled = []
            for (a, b) in tri:
                relabeled.append((a + v, b + v))
            right.append(relabeled)
        
        # Combine left and right triangulations
        for l_tri in left:
            for r_tri in right:
                # Add diagonal (0, v) and (v, n-1) if they're true diagonals
                diags = list(l_tri) + list(r_tri)
                if v > 1:  # (0, v) is a true diagonal
                    diags.append((0, v))
                if v < n - 2:  # (v, n-1) is a true diagonal
                    diags.append((v, n - 1))
                
                # Normalize: sort diagonals, ensure i < j
                normalized = []
                for (a, b) in diags:
                    if a > b:
                        a, b = b, a
                    # Only include true diagonals (not edges)
                    if b - a > 1 and not (a == 0 and b == n - 1):
                        normalized.append((a, b))
                normalized = sorted(set(normalized))
                
                if len(normalized) == n - 3:  # Correct number of diagonals
                    triangulations.append(normalized)
    
    # Remove duplicates
    unique = []
    seen = set()
    for tri in triangulations:
        key = tuple(tri)
        if key not in seen:
            seen.add(key)
            unique.append(tri)
    
    return unique


def enumerate_bcfw_cells_mhv(n, shift=(0, None)):
    """
    Enumerate BCFW cells for n-point MHV.
    
    For BCFW with shift [1, n⟩:
    - Cells correspond to BCFW diagrams with (n-3) propagators
    - For MHV, each cell is characterized by which propagator goes on-shell
    
    Args:
        n: Number of particles
        shift: Tuple (a, b) for BCFW shift [a, b⟩. Default shifts particles 0 and n-1.
    
    Returns:
        List of cells, each represented as a dict with channel information.
    """
    if n == 4:
        # For n=4, there's only one BCFW diagram
        # M_4 = M_3 × M_3 / s_12 + M_3 × M_3 / s_23
        # But for MHV, M_3 = 1, so this simplifies
        return [
            {'channel': (0, 1), 'type': 's-channel'},
            {'channel': (1, 2), 'type': 's-channel'},
        ]
    
    if n == 5:
        # For n=5, two independent channels
        return [
            {'channel': (0, 1), 'type': 's-channel'},
            {'channel': (1, 2), 'type': 's-channel'},
            {'channel': (2, 3), 'type': 's-channel'},
        ]
    
    if n == 6:
        # For n=6 MHV, we have 5 BCFW cells (Catalan(4)/2 accounting for symmetry)
        # Channels are 2-particle and 3-particle invariants
        cells = []
        
        # 2-particle channels: s_ij for adjacent pairs
        for i in range(n):
            j = (i + 1) % n
            cells.append({
                'channel': (i, j),
                'type': '2-particle',
                'particles': [i, j]
            })
        
        # 3-particle channels: s_ijk for consecutive triples
        for i in range(n):
            particles = [(i + k) % n for k in range(3)]
            cells.append({
                'channel': tuple(particles),
                'type': '3-particle',
                'particles': particles
            })
        
        # Filter to independent channels (6 + 6 = 12, but only 9 independent)
        # For BCFW [1,6⟩, relevant channels are those containing particle 1 or 6
        return cells[:9]  # Return first 9 as approximation
    
    # General case: enumerate triangulations
    tris = triangulations_of_polygon(n)
    cells = [{'triangulation': tri, 'type': 'triangulation'} for tri in tris]
    return cells


def parke_taylor_factor(tw, ordering=None):
    """
    Compute the Parke-Taylor factor for gauge theory.
    
    PT(1, 2, ..., n) = 1 / (⟨12⟩⟨23⟩...⟨n1⟩)
    
    For gravity, this appears SQUARED.
    """
    n = tw.n
    if ordering is None:
        ordering = list(range(n))
    
    result = QQ(1)
    for i in range(n):
        a = ordering[i]
        b = ordering[(i + 1) % n]
        ang = tw.angle(a, b)
        if ang == 0:
            return None
        result /= ang
    
    return result


def cell_canonical_form_mhv(cell, tw):
    """
    Compute the canonical form contribution from a single BCFW cell.
    
    For MHV gravity, each cell contributes a rational function of brackets.
    
    The form involves:
    - Parke-Taylor factors (squared for gravity)
    - Propagator factors from BCFW recursion
    - Helicity numerators
    
    This is a PLACEHOLDER implementation - the exact formula needs
    careful derivation from the BCFW recursion relations.
    """
    n = tw.n
    
    # Base: squared Parke-Taylor
    pt = parke_taylor_factor(tw)
    if pt is None:
        return None
    
    pt_squared = pt * pt
    
    # Cell-dependent factor
    cell_type = cell.get('type', '')
    
    if cell_type == '2-particle':
        # 2-particle channel: 1/s_ij factor
        i, j = cell['channel']
        s_ij = tw.mandelstam([i, j])
        if s_ij is None or s_ij == 0:
            return None
        return pt_squared / s_ij
    
    elif cell_type == '3-particle':
        # 3-particle channel: 1/s_ijk factor
        particles = cell['particles']
        s_ijk = tw.mandelstam(particles)
        if s_ijk is None or s_ijk == 0:
            return None
        return pt_squared / s_ijk
    
    else:
        # Generic cell - return Parke-Taylor squared as base
        return pt_squared


def amplitude_from_cells(tw):
    """
    Compute MHV amplitude as sum over BCFW cells.
    
    M_n^MHV = ∑_{cells} Ω_cell
    
    This should equal Hodges determinant if done correctly.
    """
    n = tw.n
    cells = enumerate_bcfw_cells_mhv(n)
    
    total = QQ(0)
    cell_contributions = []
    
    for cell in cells:
        contrib = cell_canonical_form_mhv(cell, tw)
        if contrib is None:
            return None, cells, cell_contributions
        
        cell_contributions.append({
            'cell': cell,
            'contribution': contrib
        })
        total += contrib
    
    return total, cells, cell_contributions


# =============================================================================
# GRAVITY-SPECIFIC BCFW FORMULA
# =============================================================================

def gravity_mhv_bcfw_sum(tw):
    """
    Compute 6-point MHV gravity amplitude using explicit BCFW formula.
    
    For gravity MHV, the BCFW recursion gives:
    
    M_6 = ∑_{channels} M_L × M_R × 1/P^2
    
    where M_L and M_R are lower-point amplitudes.
    For MHV, the 3-point amplitudes are simple.
    
    This is the PHYSICAL formula that should match Hodges.
    """
    n = tw.n
    if n != 6:
        raise ValueError("This formula is for n=6")
    
    # For 6-point MHV gravity, there are 9 BCFW terms
    # 6 two-particle channels + 3 three-particle channels
    
    total = QQ(0)
    
    # Get Parke-Taylor (squared for gravity)
    pt = parke_taylor_factor(tw)
    if pt is None:
        return None
    
    # Two-particle channels: s_{i,i+1}
    for i in range(n):
        ip1 = (i + 1) % n
        
        # Mandelstam s_{i,i+1}
        s = tw.mandelstam([i, ip1])
        if s is None or s == 0:
            continue
        
        # For 2-particle channel, M_3 × M_3 = 1 for MHV
        # Contribution: 1/s_{i,i+1} × (kinematic factors)
        
        # The numerator involves angle/square brackets
        # Simplified: use Parke-Taylor structure
        contrib = (pt * pt) / s
        total += contrib
    
    # Three-particle channels: s_{i,i+1,i+2}
    for i in range(n):
        ip1 = (i + 1) % n
        ip2 = (i + 2) % n
        
        s_3 = tw.mandelstam([i, ip1, ip2])
        if s_3 is None or s_3 == 0:
            continue
        
        # For 3-particle channel: M_4 × M_3 / s
        # M_3 = 1 for MHV, M_4 has specific structure
        
        # The 4-point piece involves the 4-particle amplitude
        # For MHV: M_4 = ⟨ab⟩^4/⟨12⟩⟨23⟩⟨34⟩⟨41⟩ squared for gravity
        
        # Get the 4 particles on left side
        left_particles = [i, ip1, ip2]
        right_particles = [j for j in range(n) if j not in left_particles]
        
        # Compute numerator (simplified)
        # The correct formula involves off-shell continuation
        contrib = (pt * pt) / s_3
        total += contrib
    
    return total


# =============================================================================
# TESTS
# =============================================================================

def test_triangulations():
    """Test triangulation enumeration."""
    print("Testing triangulation enumeration...")
    
    for n in range(3, 8):
        tris = triangulations_of_polygon(n)
        expected = catalan_number(n - 2)
        print(f"  n={n}: found {len(tris)} triangulations (expected {expected})")
        
        if len(tris) != expected:
            print(f"    WARNING: count mismatch!")
    
    print("Triangulation tests done!")


def test_bcfw_cells():
    """Test BCFW cell enumeration."""
    print("\nTesting BCFW cell enumeration...")
    
    for n in [4, 5, 6]:
        cells = enumerate_bcfw_cells_mhv(n)
        print(f"  n={n}: {len(cells)} cells")
        for i, cell in enumerate(cells[:3]):
            print(f"    Cell {i}: {cell}")
    
    print("BCFW cell tests done!")


def test_amplitude_structure():
    """Test amplitude computation structure."""
    print("\nTesting amplitude computation...")
    
    from src.amplituhedron.momentum_twistor import MomentumTwistorData
    
    # Try multiple seeds to find non-singular kinematics
    for seed in range(10):
        tw = MomentumTwistorData(n=6, seed=seed)
        
        if tw.is_singular():
            continue
        
        amp, cells, contribs = amplitude_from_cells(tw)
        
        if amp is None:
            continue
        
        print(f"  Seed {seed}: amplitude = {amp}")
        print(f"    Cells: {len(cells)}")
        print(f"    Non-zero contributions: {len([c for c in contribs if c['contribution'] != 0])}")
        break
    
    print("Amplitude structure tests done!")


if __name__ == '__main__':
    test_triangulations()
    test_bcfw_cells()
    test_amplitude_structure()


