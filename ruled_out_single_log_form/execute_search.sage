#!/usr/bin/env sage
# Execute search - wrapper to handle execution and monitoring
load('54.sage')

import sys
import time

def main_with_monitoring():
    """Run main with progress monitoring."""
    print("="*70)
    print("STARTING SEARCH FOR DIM-1 CANDIDATE")
    print("="*70)
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Configuration:")
    print(f"  INVARIANT_MODE: {INVARIANT_MODE}")
    print(f"  MULTI_STRATEGY_SEARCH: {MULTI_STRATEGY_SEARCH}")
    print(f"  INTERSECT_BOUNDARY_MODE: {INTERSECT_BOUNDARY_MODE}")
    print(f"  TOTAL_TRIALS: {TOTAL_TRIALS}")
    print(f"  BASE_SEED: {BASE_SEED}")
    print("="*70)
    
    try:
        main()
        print("\n" + "="*70)
        print("SEARCH COMPLETED")
        print("="*70)
    except KeyboardInterrupt:
        print("\n" + "="*70)
        print("SEARCH INTERRUPTED BY USER")
        print("="*70)
        print("Progress saved in checkpoint. Set RESUME_FROM_CHECKPOINT=True to continue.")
        sys.exit(1)
    except Exception as e:
        print("\n" + "="*70)
        print(f"ERROR: {e}")
        print("="*70)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main_with_monitoring()




