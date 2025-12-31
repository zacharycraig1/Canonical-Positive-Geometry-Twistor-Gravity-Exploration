# debug_eval.sage
import sys, os
from sage.all import *

load("src/dcp_common.sage")
load("src/kinematics_map.sage") # Explicit load
load("src/dcp_eval.sage")
load("src/spinor_sampling.sage")

def main():
    candidate_file = "dcp_cache/boundary_12_S3xS3_s45.sobj"
    Vinv_final = load(candidate_file)
    cand_vec = Vinv_final[0]
    
    data = load(cache_path("OS3_data.sobj"))
    C, triples, Vbasis = data['C'], data['triples'], data['Vbasis']
    
    # Print first few non-zero coeffs
    print("Candidate Non-Zeros:")
    nz = [(i, c) for i, c in cand_vec.dict().items()]
    for idx, coeff in nz[:10]:
        print(f"  idx={idx}, coeff={coeff}, triple={triples[idx]}")
        ti, tj, tk = triples[idx]
        print(f"    C[{ti}]={C[ti]}, C[{tj}]={C[tj]}, C[{tk}]={C[tk]}")

    # Sample
    sample = sample_spinor_helicity_conserving(n=6, seed=42)
    lambdas, tilde_lambdas = sample
    
    # Manual Trace of eval
    print("\nTracing Evaluation:")
    channels_0 = spinors_to_channels(lambdas, tilde_lambdas)
    
    # Check conversion
    channels_1 = {}
    n = 6
    def to_canon_1(key_0):
        s_1 = set([x + 1 for x in key_0])
        if len(s_1) == 3:
            comp = set(range(1, n+1)) - s_1
            return min(tuple(sorted(s_1)), tuple(sorted(comp)))
        else:
            return tuple(sorted(s_1))

    for k0, val in channels_0.items():
        k1 = to_canon_1(k0)
        channels_1[k1] = val
        
    # Check if C entries are in channels_1
    missing = []
    for ch in C:
        if ch not in channels_1:
            missing.append(ch)
    print(f"Missing channels in dict: {len(missing)}")
    if missing:
        print(f"  First missing: {missing[0]}")

    # Eval sum
    total = QQ(0)
    for idx, coeff in nz:
        ti, tj, tk = triples[idx]
        ch_i, ch_j, ch_k = C[ti], C[tj], C[tk]
        si = channels_1.get(ch_i)
        sj = channels_1.get(ch_j)
        sk = channels_1.get(ch_k)
        
        term = coeff / (si * sj * sk)
        print(f"  Term {idx}: {coeff} / ({si:.4f} * {sj:.4f} * {sk:.4f}) = {term:.4f}")
        total += term
        
    print(f"Total: {total}")

if __name__ == "__main__":
    main()

