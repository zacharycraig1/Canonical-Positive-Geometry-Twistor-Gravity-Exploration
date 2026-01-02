# run_step_0_precompute.sage
import sys, os, time
from sage.all import *

# Load modules
load("src/dcp_common.sage")
load("src/dcp_flats.sage")

def main():
    print("=" * 70)
    print("STEP 0: PRECOMPUTE FLATS AND INCOMPATIBILITY")
    print("=" * 70)
    
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)
        
    # Build matroid
    log("Building M6 matroid...")
    C, M6 = build_M6_matroid()
    m = len(C)
    log(f"  {m} channels, rank {M6.rank()}")
    
    # Build OS3
    log("Building OS3 space...")
    triples, Vbasis = build_OS3_data(C, M6)
    log(f"  OS3 dim = {Vbasis.nrows()}")
    
    # Save OS3 data for later steps
    save({'C': C, 'triples': triples, 'Vbasis': Vbasis}, cache_path("OS3_data.sobj"))
    log("Saved OS3 data.")

    # Build/Load Flats
    log("Building/Loading connected flats...")
    Gm, full_mask = get_connected_flats(M6)
    Gset = set(int(x) for x in Gm)

    # Build/Load Incompatibility
    log("Building/Loading incompatibility...")
    incompatible = get_incompatible_pairs(M6, Gm, Gset)
    
    # Flats/Incompat are automatically cached by dcp_flats functions
    log("Precomputation complete.")

if __name__ == "__main__":
    main()










