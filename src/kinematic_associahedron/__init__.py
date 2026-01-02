# src/kinematic_associahedron/__init__.py
"""
Kinematic Associahedron Module
==============================

This module implements the positive geometry for scattering amplitudes
in kinematic space, following Arkani-Hamed et al. (arXiv:1711.09102).

The key insight is that positive geometry for amplitudes lives in
KINEMATIC SPACE (Mandelstam invariants), not moduli space.

For n particles:
- Bi-adjoint scalar: Associahedron A_{n-3}
- Yang-Mills: Same geometry, different canonical form
- Gravity: KLT double copy = (Associahedron × Associahedron) / kernel

Submodules:
- associahedron: The kinematic associahedron polytope
- canonical_form: Canonical form computation
- klt_geometry: KLT kernel as geometric constraint
- gravity_from_klt: Gravity positive geometry via double copy
"""

