#!/usr/bin/env sage
"""
Quick Start: Run BCFW = Hodges Verification
============================================

This is the main entry point for Phase 2 verification.

Usage:
    sage src/amplituhedron/run_verification.sage
"""

from sage.all import *
import sys
import os

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

print("="*70)
print("GRAVITY AMPLITUHEDRON VERIFICATION")
print("="*70)
print()

# Step 1: Test momentum twistor infrastructure
print("Step 1: Testing Momentum Twistor Infrastructure...")
print("-"*50)

try:
    from src.amplituhedron.momentum_twistor import test_basic, test_positivity_search
    test_basic()
    test_positivity_search()
    print("[OK] Momentum twistor infrastructure working")
except Exception as e:
    print(f"[ERROR] Momentum twistor: {e}")
    sys.exit(1)

# Step 2: Test BCFW cells
print("\nStep 2: Testing BCFW Cell Enumeration...")
print("-"*50)

try:
    from src.amplituhedron.bcfw_cells import test_triangulations, test_bcfw_cells
    test_triangulations()
    test_bcfw_cells()
    print("[OK] BCFW cell enumeration working")
except Exception as e:
    print(f"[ERROR] BCFW cells: {e}")
    sys.exit(1)

# Step 3: Run main verification
print("\nStep 3: Running BCFW = Hodges Verification...")
print("-"*50)

try:
    from src.amplituhedron.verify_bcfw_hodges import main
    main()
except Exception as e:
    print(f"[ERROR] Verification: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "="*70)
print("VERIFICATION COMPLETE")
print("="*70)


