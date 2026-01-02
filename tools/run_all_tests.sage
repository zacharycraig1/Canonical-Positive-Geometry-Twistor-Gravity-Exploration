import sys
import os
import subprocess

def run_test(script_path):
    print(f"Running {script_path}...")
    try:
        # Use subprocess to run sage script
        # Assuming 'sage' is in path, or we use sys.executable if it's sage's python
        # Better to just use os.system for simplicity in this context, or subprocess
        ret = os.system(f"sage {script_path}")
        if ret != 0:
            print(f"FAILED: {script_path}")
            return False
        else:
            print(f"PASSED: {script_path}")
            return True
    except Exception as e:
        print(f"ERROR running {script_path}: {e}")
        return False

def main():
    scripts = [
        "tests/test_channel_identities.sage",
        "tests/test_hodges_invariance.sage",
        "tests/test_klt_equals_hodges.sage"
    ]
    
    all_passed = True
    for s in scripts:
        if not run_test(s):
            all_passed = False
            
    if all_passed:
        print("\nALL TESTS PASSED.")
        sys.exit(0)
    else:
        print("\nSOME TESTS FAILED.")
        sys.exit(1)

if __name__ == "__main__":
    main()










