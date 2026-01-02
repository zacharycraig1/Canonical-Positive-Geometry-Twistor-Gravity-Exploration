#!/usr/bin/env sage
"""
Positive Twistor Search: Find configurations in Gr_+(4,6)
==========================================================

The positive Grassmannian Gr_+(k,n) is the region where all ordered minors
are positive. For momentum twistors (k=4, n=6), we need:

    <i i+1 j j+1> > 0 for all cyclically ordered indices

Random search fails because Gr_+ is measure zero in the full Grassmannian.
Instead, we use structured generation:

1. Start with a known positive configuration
2. Use positroid cell decomposition
3. Or generate from positive parameterization

For the amplituhedron, positive configurations arise naturally from
certain BCFW-like decompositions.
"""

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


def generate_positive_twistors_structured(n=6, t_values=None):
    """
    Generate positive momentum twistors using structured parameterization.
    
    Method: Use the "positive parameterization" where twistors are
    generated from a sequence of positive parameters.
    
    For n=6, we parameterize Z_i = (1, a_i, b_i, c_i) with increasing
    parameters chosen to ensure all minors are positive.
    
    Args:
        n: Number of particles (default 6)
        t_values: Optional list of positive parameters
    
    Returns:
        List of momentum twistors guaranteed to be positive
    """
    if t_values is None:
        # Default: exponentially growing positive parameters
        t_values = [QQ(2**i) for i in range(n)]
    
    # Simple positive construction: Vandermonde-like
    # Z_i = (1, t_i, t_i^2, t_i^3)
    # This gives positive minors because Vandermonde determinants are positive
    # for distinct positive t_i
    
    Z = []
    for i in range(n):
        t = t_values[i]
        z = vector(QQ, [1, t, t**2, t**3])
        Z.append(z)
    
    return Z


def generate_positive_twistors_random_positive():
    """
    Generate positive twistors using only positive components.
    
    If all components are positive, the minors may still be positive
    with careful ordering.
    """
    n = 6
    
    # Generate with strictly positive integer components, carefully ordered
    Z = []
    for i in range(n):
        # Each successive twistor should be "larger" in some sense
        base = 10 ** i
        z = vector(QQ, [
            QQ(base + randint(1, 5)),
            QQ(base * 2 + randint(1, 5)),
            QQ(base * 3 + randint(1, 5)),
            QQ(base * 4 + randint(1, 5))
        ])
        Z.append(z)
    
    return Z


def check_positivity(Z, verbose=False):
    """
    Check if momentum twistors are in positive Grassmannian.
    
    All ordered 4-brackets <i i+1 j j+1> must be positive.
    """
    n = len(Z)
    all_positive = True
    negative_brackets = []
    
    for i in range(n):
        ip1 = (i + 1) % n
        for j in range(i + 2, n):
            jp1 = (j + 1) % n
            if jp1 == i:
                continue
            
            # 4-bracket with proper sign
            indices = [i, ip1, j, jp1]
            sorted_indices = sorted(indices)
            M = matrix([Z[k] for k in sorted_indices])
            base_det = M.det()
            
            # Sign from permutation
            inversions = sum(1 for a in range(4) for b in range(a+1, 4) 
                           if indices[a] > indices[b])
            sign = (-1) ** inversions
            signed_det = sign * base_det
            
            if signed_det <= 0:
                all_positive = False
                negative_brackets.append((i, ip1, j, jp1, signed_det))
                if verbose:
                    print(f"  Negative: <{i} {ip1} {j} {jp1}> = {signed_det}")
    
    return all_positive, negative_brackets


def twistors_to_spinors(Z):
    """Convert momentum twistors to spinor-helicity variables."""
    n = len(Z)
    
    lambdas = [vector(QQ, [z[0], z[1]]) for z in Z]
    
    def angle(i, j):
        return Z[i][0] * Z[j][1] - Z[i][1] * Z[j][0]
    
    tilde_lambdas = []
    for i in range(n):
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        
        mu_im1 = vector(QQ, [Z[im1][2], Z[im1][3]])
        mu_i = vector(QQ, [Z[i][2], Z[i][3]])
        mu_ip1 = vector(QQ, [Z[ip1][2], Z[ip1][3]])
        
        ang_i_ip1 = angle(i, ip1)
        ang_ip1_im1 = angle(ip1, im1)
        ang_im1_i = angle(im1, i)
        
        denom = ang_im1_i * ang_i_ip1
        if denom == 0:
            tilde_lambdas.append(None)
        else:
            num = mu_im1 * ang_i_ip1 + mu_i * ang_ip1_im1 + mu_ip1 * ang_im1_i
            tilde_lambdas.append(num / denom)
    
    return lambdas, tilde_lambdas


def compute_amplitude(Z):
    """Compute Hodges amplitude from momentum twistors."""
    lambdas, tilde_lambdas = twistors_to_spinors(Z)
    
    if any(t is None for t in tilde_lambdas):
        return None, "singular_tilde_lambda"
    
    return hodges_npt_mhv_canonical(lambdas, tilde_lambdas, (0, 1))


def search_positive_with_perturbation(base_Z, n_attempts=1000, epsilon=QQ(1)/QQ(10)):
    """
    Search for positive configurations by perturbing a base configuration.
    """
    print(f"Searching with perturbation from base (epsilon={epsilon})...")
    
    positive_found = []
    
    for attempt in range(n_attempts):
        # Perturb base configuration
        Z_perturbed = []
        for z in base_Z:
            z_new = vector(QQ, [
                z[i] + QQ(randint(-5, 5)) * epsilon
                for i in range(4)
            ])
            Z_perturbed.append(z_new)
        
        is_positive, _ = check_positivity(Z_perturbed)
        
        if is_positive:
            positive_found.append(Z_perturbed)
            if len(positive_found) >= 10:
                break
    
    print(f"  Found {len(positive_found)} positive configurations")
    return positive_found


def search_positive_with_scaling(n_attempts=1000):
    """
    Search for positive configurations using scaled parameters.
    
    Try different parameter scalings for the Vandermonde-like construction.
    """
    print("Searching with parameter scaling...")
    
    positive_found = []
    
    for attempt in range(n_attempts):
        set_random_seed(attempt)
        
        # Random positive base
        base = QQ(randint(2, 10))
        
        # Slightly perturbed exponential growth
        t_values = []
        t = QQ(1)
        for i in range(6):
            t_values.append(t)
            factor = base + QQ(randint(-1, 3)) / QQ(10)
            if factor <= 0:
                factor = base
            t *= factor
        
        Z = generate_positive_twistors_structured(6, t_values)
        
        is_positive, _ = check_positivity(Z)
        
        if is_positive:
            positive_found.append((t_values, Z))
            print(f"  Found positive at attempt {attempt}: t = {[float(t) for t in t_values[:3]]}...")
            
            if len(positive_found) >= 10:
                break
    
    print(f"  Found {len(positive_found)} positive configurations")
    return positive_found


def analyze_positive_configuration(Z, verbose=True):
    """
    Detailed analysis of a positive configuration.
    """
    if verbose:
        print("\n" + "="*60)
        print("POSITIVE CONFIGURATION ANALYSIS")
        print("="*60)
        
        print("\nMomentum twistors:")
        for i, z in enumerate(Z):
            print(f"  Z_{i} = {z}")
    
    # Check positivity
    is_positive, neg_brackets = check_positivity(Z, verbose=True)
    if verbose:
        print(f"\nIn positive Grassmannian: {is_positive}")
        if neg_brackets:
            print(f"  Negative brackets: {len(neg_brackets)}")
    
    # Compute amplitude
    amp, status = compute_amplitude(Z)
    if verbose:
        print(f"\nHodges amplitude:")
        print(f"  Value: {amp}")
        print(f"  Status: {status}")
        
        if amp is not None:
            try:
                amp_float = float(amp)
                print(f"  Numerical: {amp_float:.6e}")
                print(f"  Sign: {'POSITIVE' if amp_float > 0 else 'NEGATIVE'}")
            except:
                pass
    
    # Check momentum conservation
    lambdas, tilde_lambdas = twistors_to_spinors(Z)
    if all(t is not None for t in tilde_lambdas):
        P = matrix(QQ, 2, 2)
        for i in range(6):
            for a in range(2):
                for b in range(2):
                    P[a, b] += lambdas[i][a] * tilde_lambdas[i][b]
        
        if verbose:
            print(f"\nMomentum conservation: {'YES' if P == 0 else 'NO'}")
    
    return {
        'Z': Z,
        'is_positive': is_positive,
        'amplitude': amp,
        'status': status
    }


def generate_cyclic_positive_twistors(n=6):
    """
    Generate positive twistors using cyclic parameterization.
    
    For cyclic positivity, we need the twistors to form a "positive loop"
    where the cyclic 4-brackets are positive.
    
    Key insight: Use totally positive matrices (all minors positive).
    For a 4×n matrix, we can use the parameterization from
    Postnikov's positroid cells.
    """
    # A simple cyclic-positive configuration uses carefully chosen parameters
    # Based on amplituhedron literature: use positive parameters in a 
    # specific pattern
    
    # Try: Chebyshev-like oscillating pattern
    Z = []
    for i in range(n):
        theta = pi * i / n
        # Use approximate trigonometric values as rationals
        c = cos(theta)
        s = sin(theta)
        # Scale up to get nice rationals
        scale = 100
        z = vector(QQ, [
            QQ(round(scale * (1 + cos(theta)))),
            QQ(round(scale * (2 + sin(theta)))),
            QQ(round(scale * (3 + cos(2*theta)))),
            QQ(round(scale * (4 + sin(2*theta))))
        ])
        Z.append(z)
    
    return Z


def generate_amplituhedron_positive_twistors(n=6):
    """
    Generate twistors that are in the positive region by construction.
    
    Method: Use the recursive construction from Arkani-Hamed & Trnka.
    Start with a positive 4×4 matrix and extend.
    """
    # For k=4, n=6, we need a 4×6 matrix with all ordered maximal minors positive
    
    # Easiest: use a totally positive matrix
    # A simple totally positive 4×6 matrix can be constructed as follows:
    # Use Pascal's triangle or binomial coefficients
    
    Z = []
    for j in range(n):
        z = vector(QQ, [
            binomial(j, 0),
            binomial(j+1, 1),
            binomial(j+2, 2),
            binomial(j+3, 3)
        ])
        Z.append(z)
    
    return Z


def search_with_binomial_perturbation(n_attempts=1000):
    """
    Search starting from binomial matrix with perturbations.
    """
    print("Searching with binomial base + perturbations...")
    
    base_Z = generate_amplituhedron_positive_twistors(6)
    print("Base (binomial):")
    is_pos, neg = check_positivity(base_Z, verbose=True)
    print(f"  Positive: {is_pos}")
    
    if is_pos:
        return [(None, base_Z)]
    
    # Try perturbations
    positive_found = []
    for attempt in range(n_attempts):
        set_random_seed(attempt)
        
        Z_pert = []
        for z in base_Z:
            z_new = vector(QQ, [
                z[i] + QQ(randint(0, 5))
                for i in range(4)
            ])
            Z_pert.append(z_new)
        
        is_pos, _ = check_positivity(Z_pert)
        
        if is_pos:
            positive_found.append((attempt, Z_pert))
            print(f"  Found positive at attempt {attempt}")
            if len(positive_found) >= 10:
                break
    
    return positive_found


def main():
    """Main search routine."""
    print("="*70)
    print("POSITIVE MOMENTUM TWISTOR SEARCH")
    print("="*70)
    
    # Method 1: Binomial (totally positive) construction
    print("\n" + "-"*50)
    print("Method 1: Binomial (totally positive) matrix")
    print("-"*50)
    
    Z_binomial = generate_amplituhedron_positive_twistors(6)
    print("Z (binomial):")
    for i, z in enumerate(Z_binomial):
        print(f"  Z_{i} = {z}")
    
    is_pos, neg = check_positivity(Z_binomial, verbose=True)
    print(f"\nBinomial positive: {is_pos}")
    
    if is_pos:
        result = analyze_positive_configuration(Z_binomial)
        if result['amplitude'] is not None:
            try:
                amp_float = float(result['amplitude'])
                print("\n" + "="*60)
                print("KEY FINDING: POSITIVE CONFIGURATION")
                print("="*60)
                print(f"Twistors in Gr_+: {result['is_positive']}")
                print(f"Amplitude: {amp_float:.6e}")
                print(f"Sign: {'POSITIVE' if amp_float > 0 else 'NEGATIVE'}")
            except:
                pass
    
    # Method 2: Chebyshev-like cyclic
    print("\n" + "-"*50)
    print("Method 2: Chebyshev-like cyclic")
    print("-"*50)
    
    Z_chebyshev = generate_cyclic_positive_twistors(6)
    print("Z (Chebyshev):")
    for i, z in enumerate(Z_chebyshev):
        print(f"  Z_{i} = {z}")
    
    is_pos, neg = check_positivity(Z_chebyshev, verbose=True)
    print(f"\nChebyshev positive: {is_pos}")
    
    if is_pos:
        analyze_positive_configuration(Z_chebyshev)
    
    # Method 3: Binomial + perturbation
    print("\n" + "-"*50)
    print("Method 3: Binomial + perturbation search")
    print("-"*50)
    
    positive_found = search_with_binomial_perturbation(2000)
    
    if positive_found:
        print(f"\nFound {len(positive_found)} positive configurations!")
        
        # Analyze first few
        for i, (seed, Z) in enumerate(positive_found[:3]):
            print(f"\n--- Configuration {i+1} (seed={seed}) ---")
            result = analyze_positive_configuration(Z)
            
            if result['amplitude'] is not None:
                try:
                    amp_float = float(result['amplitude'])
                    print(f"Amplitude sign: {'POSITIVE' if amp_float > 0 else 'NEGATIVE'}")
                except:
                    pass
    
    return positive_found


if __name__ == "__main__":
    results = main()

