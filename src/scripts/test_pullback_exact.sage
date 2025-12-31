import sys
import os

# Ensure src is in path
sys.path.append(os.getcwd())

def run_tests():
    print("Running Exact Pullback Tests (n=4, 5, 6)...")
    
    # Load the scripts
    # Note: load() is a Sage built-in. It executes the file in the current scope.
    print("Loading test modules...")
    load("src/scripts/physics_pullback_n4.sage")
    load("src/scripts/physics_pullback_n5.sage")
    load("src/scripts/physics_pullback_n6.sage")
    
    print("-" * 40)
    
    try:
        if not physics_pullback_n4_exact():
            print("n=4 check FAILED.")
            exit(1)
            
        print("-" * 40)
        
        if not physics_pullback_n5_exact():
            print("n=5 check FAILED.")
            exit(1)
            
        print("-" * 40)
        
        if not physics_pullback_n6_exact():
            print("n=6 check FAILED.")
            exit(1)
            
        print("-" * 40)
        print("ALL TESTS PASSED")
        
    except Exception as e:
        print(f"An error occurred during testing: {e}")
        import traceback
        traceback.print_exc()
        exit(1)

if __name__ == "__main__":
    run_tests()

