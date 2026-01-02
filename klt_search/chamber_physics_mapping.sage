#!/usr/bin/env sage
"""
Chamber Physics Mapping for KLT Kernel Signature

This script builds a comprehensive atlas of signature chambers and
maps each chamber to physical scattering configurations.

Key findings to verify:
- Signature (3,3) is the mode (~70% of samples)
- Signature jumps at chamber walls (specific s_ij = 0)
- Physical interpretation of each signature class
"""
from sage.all import *
import itertools
import hashlib
import sys
import os

sys.path.insert(0, os.getcwd())


def sample_spinors(seed=None):
    """Sample momentum-conserving spinor helicity variables."""
    load('src/spinor_sampling.sage')
    return sample_spinor_helicity_conserving(n=6, seed=seed)


def ang_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def sq_bracket(tilde_lambdas, i, j):
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def mandelstam(lambdas, tilde_lambdas, i, j):
    return ang_bracket(lambdas, i, j) * sq_bracket(tilde_lambdas, i, j)


def compute_all_mandelstams(lambdas, tilde_lambdas):
    """Compute all 2-particle and 3-particle Mandelstam invariants."""
    invariants = {}
    
    # 2-particle: s_ij for i < j
    for i in range(6):
        for j in range(i+1, 6):
            invariants[(i, j)] = mandelstam(lambdas, tilde_lambdas, i, j)
    
    # 3-particle: s_ijk (cyclic)
    for i in range(6):
        j = (i + 1) % 6
        k = (i + 2) % 6
        s_ij = invariants.get((min(i,j), max(i,j)), mandelstam(lambdas, tilde_lambdas, i, j))
        s_jk = invariants.get((min(j,k), max(j,k)), mandelstam(lambdas, tilde_lambdas, j, k))
        s_ik = invariants.get((min(i,k), max(i,k)), mandelstam(lambdas, tilde_lambdas, i, k))
        invariants[(i, j, k)] = s_ij + s_jk + s_ik
    
    return invariants


def compute_klt_kernel(lambdas, tilde_lambdas):
    """Compute the 6x6 KLT kernel matrix."""
    load('src/klt.sage')
    
    class Adapter:
        def __init__(self, lam, til):
            self.lambdas = lam
            self.tilde_lambdas = til
        def get_angle(self, i, j):
            return ang_bracket(self.lambdas, i, j)
        def get_square(self, i, j):
            return sq_bracket(self.tilde_lambdas, i, j)
    
    adapter = Adapter(lambdas, tilde_lambdas)
    
    def mandelstam_func(tw, i, j):
        return tw.get_angle(i, j) * tw.get_square(i, j)
    
    perms = sorted(list(itertools.permutations([1, 2, 3])))
    
    S = matrix(QQ, 6, 6)
    for i, alpha in enumerate(perms):
        for j, beta in enumerate(perms):
            val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_func)
            S[i, j] = val if val is not None else QQ(0)
    
    return S


def get_signature(S):
    """Compute the signature (n_pos, n_neg, n_zero) of a symmetric matrix."""
    S_sym = (S + S.transpose()) / 2
    try:
        eigs = S_sym.change_ring(RDF).eigenvalues()
        n_pos = sum(1 for e in eigs if e > 1e-10)
        n_neg = sum(1 for e in eigs if e < -1e-10)
        n_zero = sum(1 for e in eigs if abs(e) <= 1e-10)
        return (n_pos, n_neg, n_zero)
    except:
        return None


def get_chamber_id(invariants):
    """
    Compute a chamber ID based on the signs of Mandelstam invariants.
    
    Chamber = region where all s_ij have definite sign.
    """
    signs = []
    for key in sorted(invariants.keys()):
        val = invariants[key]
        if val > 0:
            signs.append('+')
        elif val < 0:
            signs.append('-')
        else:
            signs.append('0')
    
    sign_str = ''.join(signs)
    chamber_hash = hashlib.md5(sign_str.encode()).hexdigest()[:8]
    
    return chamber_hash, tuple(signs)


def classify_physical_region(invariants, signs):
    """
    Classify the physical region based on Mandelstam signs.
    
    Physical scattering types:
    - Euclidean: All s_ij < 0
    - Physical 2->4: One s_ij > 0 (center of mass energy), others < 0
    - Mixed: Some positive, some negative
    """
    # Get 2-particle invariants (first 15)
    two_particle = []
    for i in range(6):
        for j in range(i+1, 6):
            two_particle.append(invariants[(i, j)])
    
    n_pos_2 = sum(1 for s in two_particle if s > 0)
    n_neg_2 = sum(1 for s in two_particle if s < 0)
    
    if n_neg_2 == 15:
        return "Euclidean"
    elif n_pos_2 == 1:
        # Find which pair is positive
        for i in range(6):
            for j in range(i+1, 6):
                if invariants[(i, j)] > 0:
                    return f"2->4 incoming ({i},{j})"
    elif n_pos_2 == 2:
        return "Physical (2 channels)"
    elif n_pos_2 == 3:
        return "Physical (3 channels)"
    else:
        return f"Mixed ({n_pos_2}+ / {n_neg_2}-)"


def build_chamber_atlas(num_samples=10000):
    """
    Build a comprehensive atlas of signature chambers.
    
    For each chamber:
    - Record signature
    - Physical classification
    - Count of samples
    """
    print("=" * 70)
    print(f"BUILDING CHAMBER ATLAS ({num_samples} samples)")
    print("=" * 70)
    
    chambers = {}  # chamber_id -> {signature, count, physical_class, signs}
    signature_counts = {}  # signature -> count
    physical_counts = {}  # physical_class -> count
    
    valid_samples = 0
    
    for sample in range(num_samples):
        if sample % 1000 == 0:
            print(f"Progress: {sample}/{num_samples}...")
        
        result = sample_spinors(seed=sample)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        # Compute invariants
        invariants = compute_all_mandelstams(lambdas, tilde_lambdas)
        
        # Check for singularities (any invariant = 0)
        if any(v == 0 for v in invariants.values()):
            continue
        
        # Compute KLT kernel and signature
        S = compute_klt_kernel(lambdas, tilde_lambdas)
        sig = get_signature(S)
        
        if sig is None:
            continue
        
        valid_samples += 1
        
        # Get chamber ID and signs
        chamber_id, signs = get_chamber_id(invariants)
        
        # Physical classification
        phys_class = classify_physical_region(invariants, signs)
        
        # Update chamber data
        if chamber_id not in chambers:
            chambers[chamber_id] = {
                'signature': sig,
                'count': 0,
                'physical_class': phys_class,
                'signs': signs
            }
        chambers[chamber_id]['count'] += 1
        
        # Track signatures
        sig_key = sig
        signature_counts[sig_key] = signature_counts.get(sig_key, 0) + 1
        
        # Track physical classifications
        physical_counts[phys_class] = physical_counts.get(phys_class, 0) + 1
    
    # Print results
    print(f"\nValid samples: {valid_samples}")
    print(f"Unique chambers: {len(chambers)}")
    
    print("\n" + "=" * 70)
    print("SIGNATURE DISTRIBUTION")
    print("=" * 70)
    
    for sig, count in sorted(signature_counts.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / valid_samples
        print(f"  {sig}: {count} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("PHYSICAL REGION DISTRIBUTION")
    print("=" * 70)
    
    for phys, count in sorted(physical_counts.items(), key=lambda x: -x[1])[:20]:
        pct = 100.0 * count / valid_samples
        print(f"  {phys}: {count} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("TOP CHAMBERS BY FREQUENCY")
    print("=" * 70)
    
    top_chambers = sorted(chambers.items(), key=lambda x: -x[1]['count'])[:20]
    for chamber_id, data in top_chambers:
        pct = 100.0 * data['count'] / valid_samples
        print(f"  Chamber {chamber_id}: {data['count']} ({pct:.1f}%)")
        print(f"    Signature: {data['signature']}")
        print(f"    Physical: {data['physical_class']}")
    
    # Analyze signature by physical region
    print("\n" + "=" * 70)
    print("SIGNATURE BY PHYSICAL REGION")
    print("=" * 70)
    
    region_signatures = {}
    for chamber_id, data in chambers.items():
        phys = data['physical_class']
        sig = data['signature']
        count = data['count']
        
        if phys not in region_signatures:
            region_signatures[phys] = {}
        region_signatures[phys][sig] = region_signatures[phys].get(sig, 0) + count
    
    for phys in sorted(region_signatures.keys()):
        print(f"\n{phys}:")
        total = sum(region_signatures[phys].values())
        for sig, count in sorted(region_signatures[phys].items(), key=lambda x: -x[1]):
            pct = 100.0 * count / total
            print(f"    {sig}: {count} ({pct:.1f}%)")
    
    # Save results
    results = {
        'chambers': chambers,
        'signature_counts': signature_counts,
        'physical_counts': physical_counts,
        'region_signatures': region_signatures,
        'valid_samples': valid_samples
    }
    
    save(results, 'klt_search/chamber_atlas_results.sobj')
    print(f"\nResults saved to klt_search/chamber_atlas_results.sobj")
    
    return results


def analyze_signature_chamber_correlation():
    """
    Analyze if specific signatures correlate with specific physical regions.
    """
    try:
        results = load('klt_search/chamber_atlas_results.sobj')
    except:
        print("Run build_chamber_atlas first!")
        return
    
    print("\n" + "=" * 70)
    print("SIGNATURE-CHAMBER CORRELATION ANALYSIS")
    print("=" * 70)
    
    region_sigs = results['region_signatures']
    
    # For each signature, find the most common physical region
    signature_to_region = {}
    
    for phys, sigs in region_sigs.items():
        for sig, count in sigs.items():
            if sig not in signature_to_region:
                signature_to_region[sig] = {}
            signature_to_region[sig][phys] = signature_to_region[sig].get(phys, 0) + count
    
    print("\nFor each signature, most common physical region:")
    for sig in sorted(signature_to_region.keys()):
        regions = signature_to_region[sig]
        total = sum(regions.values())
        top_region = max(regions.items(), key=lambda x: x[1])
        pct = 100.0 * top_region[1] / total
        print(f"  {sig}: {top_region[0]} ({pct:.1f}%)")


if __name__ == "__main__":
    import sys
    
    # Default to 500 samples for quick testing, use command line arg for more
    num_samples = 500
    if len(sys.argv) > 1:
        try:
            num_samples = int(sys.argv[1])
        except:
            pass
    
    print(f"Running with {num_samples} samples...")
    build_chamber_atlas(num_samples=num_samples)
    
    # Analyze correlations
    analyze_signature_chamber_correlation()

