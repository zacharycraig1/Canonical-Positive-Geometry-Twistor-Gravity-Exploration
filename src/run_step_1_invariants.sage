# run_step_1_invariants.sage
import sys, os, time
from sage.all import *

# Load modules
load("src/dcp_common.sage")
load("src/dcp_flats.sage")
load("src/dcp_invariants.sage")

def main():
    print("=" * 70)
    print("STEP 1: COMPUTE INVARIANTS")
    print("=" * 70)

    # Load OS3 data
    try:
        data = load(cache_path("OS3_data.sobj"))
        C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']
        log("Loaded OS3 data.")
    except Exception as e:
        log(f"Error loading OS3 data: {e}. Run step 0 first.")
        return

    # S3xS3
    log("Computing S3xS3 invariants...")
    V_s3xs3 = compute_boundary_stabilizer_invariants(C, triples, Vbasis, mode='S3xS3')
    save(V_s3xs3, cache_path("basis_S3xS3.sobj"))
    log(f"Saved {len(V_s3xs3)} S3xS3 invariants to basis_S3xS3.sobj")

    # S3xS3Z2
    log("Computing S3xS3Z2 invariants...")
    V_s3xs3z2 = compute_boundary_stabilizer_invariants(C, triples, Vbasis, mode='S3xS3Z2')
    save(V_s3xs3z2, cache_path("basis_S3xS3Z2.sobj"))
    log(f"Saved {len(V_s3xs3z2)} S3xS3Z2 invariants to basis_S3xS3Z2.sobj")

    # S6
    log("Computing S6 invariants...")
    V_s6 = compute_S6_invariants(C, triples, Vbasis)
    save(V_s6, cache_path("basis_S6.sobj"))
    log(f"Saved {len(V_s6)} S6 invariants to basis_S6.sobj")

if __name__ == "__main__":
    main()









