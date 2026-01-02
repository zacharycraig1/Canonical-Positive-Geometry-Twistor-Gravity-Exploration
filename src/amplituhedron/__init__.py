# src/amplituhedron/__init__.py
"""
Gravity Amplituhedron Module
============================

This module implements the BCFW cell construction for MHV gravity amplitudes.

Key components:
- momentum_twistor.sage: Momentum twistor kinematics with positivity checks
- bcfw_cells.sage: BCFW cell enumeration and canonical forms
- hodges_twistor.sage: Independent Hodges determinant implementation
- verify_bcfw_hodges.sage: Verification tests

The goal is to prove:
    BCFW cell sum = Hodges determinant

for 6-point MHV gravity, establishing a positive geometry for gravity.
"""

__version__ = "0.1.0"
__author__ = "Gravity Amplituhedron Project"


