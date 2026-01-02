#!/usr/bin/env sage
# =============================================================================
# Main Runner: Find Positive Geometry for 6-Point MHV Gravity
# =============================================================================
# Executes the complete geometry search as specified in the plan.

from sage.all import *
import sys
import os
import json
from datetime import datetime

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load the geometry finder
load("src/gravity_proof/geometry_finder.sage")


def main(seed=42, save_results=True):
    """
    Main execution: Find the positive geometry R6.
    
    Args:
        seed: Random seed for kinematics
        save_results: Whether to save results to file
        
    Returns:
        (finder, results) tuple
    """
    print("="*80)
    print("FINDING POSITIVE GEOMETRY FOR 6-POINT MHV GRAVITY")
    print("="*80)
    print(f"Date: {datetime.now().isoformat()}")
    print(f"Seed: {seed}")
    print("="*80)
    
    # Run the analysis
    finder, results = GeometryFinder(seed=seed), None
    results = finder.run_full_analysis()
    
    # Save results if requested
    if save_results:
        output_file = "geometry_search_results.json"
        
        # Convert to JSON-serializable format
        def make_serializable(obj):
            if isinstance(obj, (int, float, str, bool, type(None))):
                return obj
            elif isinstance(obj, dict):
                return {k: make_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple)):
                return [make_serializable(item) for item in obj]
            else:
                return str(obj)
        
        serializable_results = make_serializable(results)
        
        try:
            with open(output_file, 'w') as f:
                json.dump({
                    'timestamp': datetime.now().isoformat(),
                    'seed': seed,
                    'results': serializable_results
                }, f, indent=2)
            print(f"\nResults saved to: {output_file}")
        except Exception as e:
            print(f"\nWarning: Could not save results: {e}")
    
    # Print final conclusion
    print(f"\n" + "="*80)
    print("CONCLUSION")
    print("="*80)
    
    correct_option = results.get('correct_option')
    
    if correct_option:
        print(f"\n✓ Identified geometry: Option {correct_option}")
        print(f"  This option defines R6 with the correct vertex structure.")
        
        if 'symbolic_results' in results:
            sr = results['symbolic_results']
            if 'numerical_comparison' in sr:
                nc = sr['numerical_comparison']
                if nc.get('status') == 'success' and nc.get('num_matches', 0) > 0:
                    print(f"\n✓ Canonical form matches Hodges amplitude at {nc['num_matches']} solution(s)")
                    print(f"  The positive geometry approach successfully reproduces the gravity amplitude!")
                else:
                    print(f"\n⚠ Canonical form does not match Hodges numerically")
                    print(f"  Further investigation needed on the geometric derivation")
    else:
        print(f"\n⚠ Could not definitively identify correct geometry")
        print(f"  Options tested: B1, B2, B3, B4")
        print(f"  None gave the expected 6-vertex structure")
        print(f"\nPossible next steps:")
        print(f"  - Refine numerical tolerances")
        print(f"  - Test combinations of conditions")
        print(f"  - Check if additional constraints are needed")
    
    return finder, results


if __name__ == "__main__":
    # Allow command-line seed specification
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    
    finder, results = main(seed=seed)

