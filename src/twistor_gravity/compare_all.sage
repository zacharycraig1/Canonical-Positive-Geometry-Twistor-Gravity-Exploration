#!/usr/bin/env sage
"""
Cross-Validation of All Approaches for MHV Gravity
===================================================

This module computes the 6-point MHV gravity amplitude using all available
methods and verifies they give the same answer.

Methods compared:
1. Hodges determinant (spinor variables)
2. Hodges determinant (twistor variables)
3. CHY sum over scattering equation solutions
4. Gravituhedron/double-copy approach
5. Twistor string MHV vertex

The goal is to identify which formulation:
- Makes positivity manifest
- Gives a clean canonical form = amplitude
- Connects to known structures

Additionally, we compare the positive regions:
| Approach        | Space         | Positivity Condition            |
|-----------------|---------------|--------------------------------|
| CHY             | M_{0,6}       | z_4 > z_5 > z_6 > 0 (failed)   |
| Momentum Twistor| Gr_+(4,6)     | All <i i+1 j j+1> > 0          |
| Twistor String  | Line in PT    | Degree-1 map positivity        |
| Celestial       | Celestial S²  | OPE coefficient signs          |
"""

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))


def generate_test_kinematics(seed=42):
    """
    Generate consistent kinematics for all methods.
    
    Returns momentum twistors that can be converted to spinors.
    """
    set_random_seed(seed)
    n = 6
    
    # Generate random momentum twistors
    Z = []
    for i in range(n):
        z = vector(QQ, [QQ(randint(1, 20)) for _ in range(4)])
        Z.append(z)
    
    return Z


def twistors_to_spinors(Z):
    """
    Convert momentum twistors to spinor-helicity variables.
    
    λ_i = (Z_i^0, Z_i^1)
    λ̃_i reconstructed from incidence relations
    """
    n = len(Z)
    
    # Lambda from first two components
    lambdas = [vector(QQ, [z[0], z[1]]) for z in Z]
    
    # Angle brackets for reconstruction
    def angle(i, j):
        return Z[i][0] * Z[j][1] - Z[i][1] * Z[j][0]
    
    # Reconstruct tilde_lambda
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


def compute_hodges_spinor(lambdas, tilde_lambdas, negative_helicity=(0, 1)):
    """Compute Hodges amplitude using spinor-based method."""
    try:
        from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
        
        if any(t is None for t in tilde_lambdas):
            return None, "singular_kinematics"
        
        result, status = hodges_npt_mhv_canonical(
            lambdas, tilde_lambdas, negative_helicity
        )
        return result, status
    except Exception as e:
        return None, str(e)


def compute_hodges_twistor(Z, negative_helicity=(0, 1)):
    """Compute Hodges amplitude using twistor-based method."""
    try:
        load("src/twistor_gravity/hodges_twistor.sage")
        ht = HodgesTwistor(Z, negative_helicity=negative_helicity)
        return ht.mhv_amplitude()
    except Exception as e:
        return None, str(e)


def compute_chy_amplitude(lambdas, tilde_lambdas):
    """Compute amplitude via CHY formula with scattering equations."""
    try:
        from src.kinematics.spinors import SpinorKinematics
        
        if any(t is None for t in tilde_lambdas):
            return None, "singular_kinematics"
        
        kin = SpinorKinematics(6, lambdas, tilde_lambdas)
        
        load("src/gravity_proof/scattering_solver.sage")
        load("src/gravity_proof/psi_matrix.sage")
        
        solver = ScatteringEquationSolver(kin)
        solutions = solver.solve_numerical()
        
        if len(solutions) != 6:
            return None, f"got_{len(solutions)}_solutions"
        
        # Helicity factor
        ang_12 = kin.angle(0, 1)
        helicity_factor = ang_12**8
        
        # Sum over solutions
        chy_sum = 0
        for sol in solutions:
            z4, z5, z6 = sol['z4'], sol['z5'], sol['z6']
            
            psi = PsiMatrixMHV(kin, z4, z5, z6)
            try:
                pf_prime = psi.reduced_pfaffian_standard((0, 1))
                if pf_prime == 0:
                    pf_prime = psi.reduced_pfaffian_delete_3_6()
            except:
                try:
                    pf_prime = psi.reduced_pfaffian_delete_3_6()
                except:
                    continue
            
            J = solver.jacobian_matrix(z4, z5, z6)
            det_J = J.det()
            
            if det_J == 0:
                continue
            
            chy_sum += (pf_prime**2) / det_J
        
        return helicity_factor * chy_sum, "ok"
    except Exception as e:
        return None, str(e)


def compute_gravituhedron(Z):
    """Compute amplitude via gravituhedron/double-copy approach."""
    try:
        load("src/twistor_gravity/gravituhedron.sage")
        gc = GravituhedronCandidate(Z)
        return gc.compute_candidate_amplitude('hodges')
    except Exception as e:
        return None, str(e)


def compute_twistor_string(Z):
    """Compute amplitude via twistor string MHV vertex."""
    try:
        load("src/twistor_gravity/twistor_string.sage")
        ts = TwistorStringGravity(Z)
        amp = ts.scattering_amplitude('vertex')
        return amp, "ok" if amp is not None else "singular"
    except Exception as e:
        return None, str(e)


def check_positivity_regions(Z, lambdas, tilde_lambdas):
    """
    Check positivity conditions in each formulation.
    
    Returns dict with positivity status for each approach.
    """
    results = {}
    
    # 1. Twistor positivity: all <i i+1 j j+1> > 0
    try:
        load("src/twistor_gravity/hodges_twistor.sage")
        ht = HodgesTwistor(Z)
        results['twistor_positive'] = ht.is_positive()
    except:
        results['twistor_positive'] = None
    
    # 2. Gravituhedron positivity
    try:
        load("src/twistor_gravity/gravituhedron.sage")
        gc = GravituhedronCandidate(Z)
        results['gravituhedron_positive'] = gc.is_in_positive_grassmannian()
    except:
        results['gravituhedron_positive'] = None
    
    # 3. CHY/worldsheet positivity (would need scattering equation solutions)
    results['worldsheet_positive'] = 'requires_solution_analysis'
    
    # 4. Celestial positivity (theoretical)
    results['celestial_positive'] = 'theoretical_only'
    
    return results


def cross_validate_all(seed=42, verbose=True):
    """
    Compute MHV gravity amplitude via all methods and compare.
    
    Args:
        seed: Random seed for kinematics
        verbose: Print detailed output
    
    Returns:
        dict with all results
    """
    if verbose:
        print("="*70)
        print("CROSS-VALIDATION OF ALL MHV GRAVITY APPROACHES")
        print("="*70)
        print(f"Seed: {seed}")
    
    # Generate kinematics
    Z = generate_test_kinematics(seed)
    lambdas, tilde_lambdas = twistors_to_spinors(Z)
    
    if verbose:
        print(f"\nGenerated {len(Z)} momentum twistors")
        singular = sum(1 for t in tilde_lambdas if t is None)
        if singular > 0:
            print(f"WARNING: {singular} particles have singular tilde_lambda")
    
    results = {
        'seed': seed,
        'n': 6,
        'methods': {}
    }
    
    # Method 1: Hodges (spinor)
    if verbose:
        print("\n" + "-"*50)
        print("Method 1: Hodges (spinor variables)")
    
    amp1, status1 = compute_hodges_spinor(lambdas, tilde_lambdas)
    results['methods']['hodges_spinor'] = {
        'amplitude': amp1,
        'status': status1
    }
    if verbose:
        if amp1 is not None:
            try:
                print(f"  Amplitude: {float(amp1):.6e}")
            except:
                print(f"  Amplitude: {amp1}")
        else:
            print(f"  Failed: {status1}")
    
    # Method 2: Hodges (twistor)
    if verbose:
        print("\n" + "-"*50)
        print("Method 2: Hodges (twistor variables)")
    
    amp2, status2 = compute_hodges_twistor(Z)
    results['methods']['hodges_twistor'] = {
        'amplitude': amp2,
        'status': status2
    }
    if verbose:
        if amp2 is not None:
            try:
                print(f"  Amplitude: {float(amp2):.6e}")
            except:
                print(f"  Amplitude: {amp2}")
        else:
            print(f"  Failed: {status2}")
    
    # Method 3: CHY
    if verbose:
        print("\n" + "-"*50)
        print("Method 3: CHY (scattering equations)")
    
    amp3, status3 = compute_chy_amplitude(lambdas, tilde_lambdas)
    results['methods']['chy'] = {
        'amplitude': amp3,
        'status': status3
    }
    if verbose:
        if amp3 is not None:
            try:
                print(f"  Amplitude: {float(amp3):.6e}")
            except:
                print(f"  Amplitude: {amp3}")
        else:
            print(f"  Failed: {status3}")
    
    # Method 4: Gravituhedron
    if verbose:
        print("\n" + "-"*50)
        print("Method 4: Gravituhedron (double copy)")
    
    amp4, status4 = compute_gravituhedron(Z)
    results['methods']['gravituhedron'] = {
        'amplitude': amp4,
        'status': status4
    }
    if verbose:
        if amp4 is not None:
            try:
                print(f"  Amplitude: {float(amp4):.6e}")
            except:
                print(f"  Amplitude: {amp4}")
        else:
            print(f"  Failed: {status4}")
    
    # Method 5: Twistor string
    if verbose:
        print("\n" + "-"*50)
        print("Method 5: Twistor string (MHV vertex)")
    
    amp5, status5 = compute_twistor_string(Z)
    results['methods']['twistor_string'] = {
        'amplitude': amp5,
        'status': status5
    }
    if verbose:
        if amp5 is not None:
            try:
                print(f"  Amplitude: {float(amp5):.6e}")
            except:
                print(f"  Amplitude: {amp5}")
        else:
            print(f"  Failed: {status5}")
    
    # Check positivity
    if verbose:
        print("\n" + "="*70)
        print("POSITIVITY ANALYSIS")
        print("="*70)
    
    positivity = check_positivity_regions(Z, lambdas, tilde_lambdas)
    results['positivity'] = positivity
    
    if verbose:
        for region, is_pos in positivity.items():
            print(f"  {region}: {is_pos}")
    
    # Compare results
    if verbose:
        print("\n" + "="*70)
        print("COMPARISON")
        print("="*70)
    
    # Get reference value (Hodges spinor)
    ref_amp = amp1
    ref_name = 'hodges_spinor'
    
    if ref_amp is not None:
        comparisons = []
        for name, data in results['methods'].items():
            amp = data['amplitude']
            if amp is not None and name != ref_name:
                try:
                    amp_float = float(amp)
                    ref_float = float(ref_amp)
                    ratio = amp_float / ref_float if ref_float != 0 else float('inf')
                    rel_diff = abs(amp_float - ref_float) / abs(ref_float) if ref_float != 0 else float('inf')
                    match = rel_diff < 1e-8
                    
                    comparisons.append({
                        'method': name,
                        'ratio': ratio,
                        'rel_diff': rel_diff,
                        'match': match
                    })
                    
                    if verbose:
                        status = 'MATCH' if match else 'MISMATCH'
                        print(f"  {name} vs {ref_name}: ratio={ratio:.6f}, rel_diff={rel_diff:.2e} [{status}]")
                except:
                    if verbose:
                        print(f"  {name}: could not compare")
        
        results['comparisons'] = comparisons
    else:
        if verbose:
            print("  Reference amplitude (Hodges spinor) failed, cannot compare")
    
    # Summary
    if verbose:
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        
        successful = sum(1 for m, d in results['methods'].items() if d['amplitude'] is not None)
        print(f"  Methods computed successfully: {successful}/{len(results['methods'])}")
        
        if 'comparisons' in results:
            matches = sum(1 for c in results['comparisons'] if c['match'])
            print(f"  Methods matching reference: {matches}/{len(results['comparisons'])}")
    
    return results


if __name__ == "__main__":
    # Run cross-validation with multiple seeds
    seeds = [42, 123, 456]
    all_results = []
    
    for seed in seeds:
        result = cross_validate_all(seed=seed)
        all_results.append(result)
        print("\n" + "#"*70 + "\n")
    
    # Final summary
    print("="*70)
    print("FINAL SUMMARY ACROSS ALL SEEDS")
    print("="*70)
    
    for i, result in enumerate(all_results):
        seed = result['seed']
        successful = sum(1 for m, d in result['methods'].items() if d['amplitude'] is not None)
        matches = sum(1 for c in result.get('comparisons', []) if c['match'])
        total = len(result.get('comparisons', []))
        print(f"  Seed {seed}: {successful}/5 computed, {matches}/{total} matched reference")

