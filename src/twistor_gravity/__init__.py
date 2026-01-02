# src/twistor_gravity/__init__.py
"""
Twistor Gravity Module
======================

This module explores twistor-space formulations of MHV gravity amplitudes
to find the positive geometry for 6-point MHV gravity.

Key approaches:
1. Momentum Twistors + Gravituhedron: Extend amplituhedron infrastructure
2. Twistor String: Implement Skinner's worldsheet-to-twistor formula
3. Celestial Holography: Map to 2D CFT on celestial sphere

Components:
- hodges_twistor.sage: Hodges formula in pure twistor variables
- gravituhedron.sage: Candidate positive geometry for gravity
- twistor_string.sage: Skinner formula implementation
- celestial_map.sage: Celestial amplitude transform
- compare_all.sage: Cross-validation of all approaches

The goal is to find which formulation makes positivity manifest and
gives a clean canonical form = amplitude.
"""

__version__ = "0.1.0"
__author__ = "Gravity Amplituhedron Project"

