
from sage.all import *
from itertools import permutations
import numpy as np
import hashlib

load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')

# --- Helper Classes & Functions ---

class SpinorHelicityAdapter:
    def __init__(self, lambdas, tilde_lambdas):
        self.n = len(lambdas)
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas
    def get_angle(self, i, j):
        return ang_bracket(self.lambdas, i, j)
    def get_square(self, i, j):
        return sq_bracket(self.tilde_lambdas, i, j)

def mandelstam_simple(tw, i, j):
    return tw.get_angle(i, j) * tw.get_square(i, j)

def get_region_signature_robust(adapter, mandelstam_func):
    # 1. Compute Invariants & Signs
    invariants = []
    
    # 2-particle (15 for n=6)
    # Ordered pairs (i,j) with i<j
    pairs = []
    for i in range(6):
        for j in range(i+1, 6):
            pairs.append((i,j))
            val = mandelstam_func(adapter, i, j)
            invariants.append(val)
            
    # 3-particle (cyclic triples)
    # s_123, s_234, s_345, s_456, s_561, s_612
    triples = []
    for i in range(6):
        j = (i+1)%6
        k = (i+2)%6
        triples.append((i,j,k))
        # s_ijk = s_ij + s_jk + s_ki
        s_ij = mandelstam_func(adapter, i, j)
        s_jk = mandelstam_func(adapter, j, k)
        s_ki = mandelstam_func(adapter, k, i)
        val = s_ij + s_jk + s_ki
        invariants.append(val)
        
    # Check for zeros (boundary)
    if any(x == 0 for x in invariants):
        return None, None, None, True # near_singular/boundary
        
    signs = tuple([1 if x > 0 else -1 for x in invariants])
    chamber_id = hashlib.md5(str(signs).encode()).hexdigest()[:8]
    
    # 2. Compute KLT/Intersection Matrix Signature
    # Use Symmetric KLT kernel
    permuted_set = [1, 2, 3] 
    basis_perms = sorted(list(permutations(permuted_set)))
    
    matrix_data = []
    for alpha in basis_perms:
        row = []
        for beta in basis_perms:
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_func)
            row.append(val)
        matrix_data.append(row)
        
    S = matrix(QQ, matrix_data)
    S_sym = (S + S.transpose()) / 2
    
    # Diagnostics
    try:
        detS = S.det() # Not symmetric one? No, use original KLT for det check?
        # Actually KLT is symmetric.
        detS = S_sym.det()
    except: detS = 0
    
    if detS == 0:
        return None, None, None, True
        
    # Eigenvalues
    try:
        evals = S_sym.change_ring(RDF).eigenvalues()
        min_abs_eig = min(abs(e) for e in evals)
        
        n_pos = sum(1 for e in evals if e > 1e-6)
        n_neg = sum(1 for e in evals if e < -1e-6)
        n_zero = sum(1 for e in evals if abs(e) <= 1e-6)
        signature = (n_pos, n_neg, n_zero)
        
        near_singular = (min_abs_eig < 1e-10) or (abs(detS) < 1e-20)
        
    except:
        return None, None, None, True
        
    return chamber_id, signs, signature, near_singular

def get_region_tags(signs):
    # signs vector length 21: 15 pairs + 6 triples
    # Pairs indices: (0,1), (0,2)...
    # Triples indices: 15, 16...
    
    # 1. Euclidean-like: All 2-particle invariants < 0
    # indices 0..14
    s2_signs = signs[:15]
    is_euclidean = all(s < 0 for s in s2_signs)
    
    # 2. 2->4 Heuristic
    # Incoming pair (a,b) > 0, all others involving a or b < 0?
    # Scan all 15 pairs
    incoming_candidates = []
    pair_idx = 0
    for i in range(6):
        for j in range(i+1, 6):
            if signs[pair_idx] > 0: # s_ab > 0
                # Check condition: s_ak < 0 and s_bk < 0 for k != a,b
                is_2to4 = True
                for k in range(6):
                    if k == i or k == j: continue
                    # find index for (min(i,k), max(i,k))
                    idx_ak = get_pair_index(i, k)
                    idx_bk = get_pair_index(j, k)
                    
                    if signs[idx_ak] > 0 or signs[idx_bk] > 0:
                        is_2to4 = False
                        break
                
                if is_2to4:
                    incoming_candidates.append((i,j))
            pair_idx += 1
            
    tags = []
    if is_euclidean: tags.append("Euclidean")
    if incoming_candidates: tags.append(f"2->4 {incoming_candidates}")
    
    return tags

def get_pair_index(i, j):
    # Map (i,j) to 0..14
    u, v = sorted((i, j))
    # Formula: index = sum_{k=0}^{u-1} (6-1-k) + (v-u-1)
    # Or just loop
    idx = 0
    for a in range(6):
        for b in range(a+1, 6):
            if a==u and b==v: return idx
            idx += 1
    return -1

def build_atlas(num_samples=200):
    print(f"Building Signature Atlas (n=6, {num_samples} samples)...")
    
    chamber_stats = {} # id -> {count, signature, tags, signs}
    signature_stats = {} # sig -> count
    
    total = 0
    kept = 0
    near_sing = 0
    
    # Try more samples to find Euclidean
    # Euclidean requires s_ij < 0.
    # In (2,2) signature (split metric), s_ij can be negative.
    # But usually physical scattering implies some s > 0.
    # To hit Euclidean, maybe we need specific sampling?
    # Random spinors cover phase space sparsely?
    # Actually, random complex spinors -> random real s_ij?
    # sample_spinor_helicity_conserving uses random integers.
    # s_ij = <ij>[ij]. For real momenta, s_ij is real.
    # If lambda, tilde are random QQ, s_ij are random QQ.
    # Probability of all 15 s_ij < 0 is small? (1/2)^15?
    # Yes, 3e-5.
    # So we need more samples or targeted sampling.
    
    extra_samples_if_needed = 1000
    
    for k in range(num_samples):
        total += 1
        try:
            lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6)
            adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
            
            cid, signs, sig, bad = get_region_signature_robust(adapter, mandelstam_simple)
            
            if bad:
                near_sing += 1
                continue
                
            kept += 1
            
            # Tags
            tags = get_region_tags(signs)
            
            # Record
            if cid not in chamber_stats:
                chamber_stats[cid] = {
                    'count': 0,
                    'signature': sig,
                    'tags': tags,
                    'signs': signs
                }
            chamber_stats[cid]['count'] += 1
            
            # Check signature consistency for chamber
            if chamber_stats[cid]['signature'] != sig:
                print(f"WARNING: Signature jump within chamber {cid}! {chamber_stats[cid]['signature']} -> {sig}")
                
            if sig not in signature_stats:
                signature_stats[sig] = 0
            signature_stats[sig] += 1
            
        except Exception:
            continue
            
    print(f"\nSampling Complete.")
    print(f"Total: {total}, Kept: {kept}, Near-Singular/Boundary: {near_sing}")
    
    print("\nTop Chambers:")
    sorted_chambers = sorted(chamber_stats.items(), key=lambda x: x[1]['count'], reverse=True)
    
    for cid, data in sorted_chambers[:10]:
        print(f"Chamber {cid}: Count={data['count']}, Sig={data['signature']}, Tags={data['tags']}")
        
    print("\nSignature Distribution:")
    for sig, count in signature_stats.items():
        print(f"  {sig}: {count}")
        
    # Check correlation: (3,3) and Euclidean
    euclidean_sigs = []
    for cid, data in chamber_stats.items():
        if "Euclidean" in data['tags']:
            euclidean_sigs.append(data['signature'])
            
    if euclidean_sigs:
        print(f"\nEuclidean Regions Signature: {set(euclidean_sigs)}")
    else:
        print("\nNo Euclidean regions found.")

if __name__ == "__main__":
    build_atlas(200)

