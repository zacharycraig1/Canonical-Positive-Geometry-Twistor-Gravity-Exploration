#!/usr/bin/env sage
# Quick environment check
print("SageMath version:", version())
print("Python version:", sys.version)
print("Working directory:", os.getcwd())
print("CPU count:", os.cpu_count())
print("\nTesting imports...")
try:
    import numpy as np
    print("  numpy: OK")
except:
    print("  numpy: FAILED")
try:
    from sage.all import *
    print("  sage.all: OK")
except:
    print("  sage.all: FAILED")
print("\nEnvironment check complete.")









