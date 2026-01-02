#!/usr/bin/env sage
"""
Kinematic Space Geometry for MHV Gravity
=========================================

This module explores positive geometry in KINEMATIC SPACE directly,
rather than twistor space or moduli space.

Key insight from research:
- For Yang-Mills, the positive geometry is the amplituhedron in Gr_+(k,n)
- For gravity, the geometry might be in a DIFFERENT space
- The double-copy structure suggests: Kinematic space with positivity

Kinematic space for n particles:
- Coordinates: Mandelstam invariants s_{ij}, s_{ijk}, etc.
- Momentum conservation: constraints reduce dimension
- Positivity: certain regions may correspond to amplitudes

For n=6:
- Independent Mandelstams: s12, s23, s34, s45, s56, s61, s123, s234, s345 (9 variables)
- With constraints: dim = 9 - 4 = 5

The gravity amplitude may be the canonical form of a positive region
in this kinematic space!
"""

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


class KinematicSpaceGeometry:
    """
    Explores positive geometry in kinematic space for MHV gravity.
    """
    
    def __init__(self, kinematics=None, seed=42):
        """
        Initialize with 6-particle kinematics.
        
        Args:
            kinematics: Optional SpinorKinematics object
            seed: Random seed if generating kinematics
        """
        if kinematics is not None:
            self.kin = kinematics
        else:
            self.kin = SpinorKinematics.random_rational(6, seed=seed)
        
        self.n = 6
        self._compute_invariants()
    
    def _compute_invariants(self):
        """Compute all kinematic invariants."""
        n = self.n
        
        # Two-particle invariants s_ij
        self.s2 = {}
        for i in range(n):
            for j in range(i+1, n):
                self.s2[(i, j)] = self.kin.s(i, j)
        
        # Three-particle invariants s_ijk
        self.s3 = {}
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    # s_ijk = s_ij + s_jk + s_ik
                    self.s3[(i, j, k)] = (
                        self.s2[(i, j)] + 
                        self.s2[(j, k)] + 
                        self.s2[(i, k)]
                    )
    
    def get_mandelstam(self, *indices):
        """Get Mandelstam invariant for given particle indices."""
        indices = tuple(sorted(indices))
        if len(indices) == 2:
            return self.s2.get(indices, QQ(0))
        elif len(indices) == 3:
            return self.s3.get(indices, QQ(0))
        else:
            # Compute from sum
            total = QQ(0)
            for i in range(len(indices)):
                for j in range(i+1, len(indices)):
                    total += self.s2.get((indices[i], indices[j]), QQ(0))
            return total
    
    def factorization_poles(self):
        """
        List all factorization poles for 6-point amplitude.
        
        These are where s_I = 0 for some subset I.
        They correspond to boundaries of the positive geometry.
        """
        poles = []
        
        # Two-particle poles
        for (i, j), s in self.s2.items():
            poles.append({
                'type': 'two_particle',
                'particles': (i, j),
                'invariant': f's_{i}{j}',
                'value': s
            })
        
        # Three-particle poles (only 3 independent for n=6)
        # s_012 = s_345, s_123 = s_450, s_234 = s_501
        for (i, j, k), s in self.s3.items():
            poles.append({
                'type': 'three_particle',
                'particles': (i, j, k),
                'invariant': f's_{i}{j}{k}',
                'value': s
            })
        
        return poles
    
    def boundary_structure(self, verbose=True):
        """
        Analyze the boundary structure of the kinematic region.
        
        For gravity, the amplitude has poles at:
        - s_ij = 0 (two-particle channels)
        - s_ijk = 0 (three-particle channels)
        
        The positive geometry should have boundaries at these poles.
        """
        poles = self.factorization_poles()
        
        if verbose:
            print("="*60)
            print("KINEMATIC SPACE BOUNDARY STRUCTURE")
            print("="*60)
            print(f"\n{len(poles)} factorization poles:")
        
        pos_count = 0
        neg_count = 0
        zero_count = 0
        
        for pole in poles:
            val = pole['value']
            try:
                val_float = float(val)
                sign = '+' if val_float > 0 else '-' if val_float < 0 else '0'
            except:
                sign = '?'
            
            if sign == '+':
                pos_count += 1
            elif sign == '-':
                neg_count += 1
            else:
                zero_count += 1
            
            if verbose:
                print(f"  {pole['invariant']} = {float(val):.4f} [{sign}]")
        
        if verbose:
            print(f"\nPositive: {pos_count}, Negative: {neg_count}, Zero: {zero_count}")
        
        return {
            'poles': poles,
            'positive_count': pos_count,
            'negative_count': neg_count,
            'zero_count': zero_count
        }
    
    def check_physical_region(self, verbose=True):
        """
        Check if kinematics are in the physical region.
        
        Physical region (for real scattering):
        - All s_ij real
        - Certain positivity/negativity constraints
        
        For gravity, the amplitude should be well-defined in this region.
        """
        if verbose:
            print("\n" + "="*60)
            print("PHYSICAL REGION CHECK")
            print("="*60)
        
        # Check if all invariants are real (QQ is always real)
        all_real = True
        for (i, j), s in self.s2.items():
            # For QQ, values are always real
            try:
                _ = float(s)
            except:
                all_real = False
                if verbose:
                    print(f"  s_{i}{j} is complex")
        
        if verbose and all_real:
            print("  All invariants are real ✓")
        
        # Check momentum conservation (s_123 = s_456)
        s123 = self.get_mandelstam(0, 1, 2)
        s456 = self.get_mandelstam(3, 4, 5)
        
        if verbose:
            print(f"\nMomentum conservation:")
            print(f"  s_123 = {float(s123):.6f}")
            print(f"  s_456 = {float(s456):.6f}")
            print(f"  Match: {s123 == s456}")
        
        return all_real and (s123 == s456)
    
    def compute_amplitude(self, verbose=True):
        """
        Compute the 6-point MHV gravity amplitude.
        
        Uses the Hodges formula as reference.
        """
        lambdas = self.kin.lambdas
        tilde_lambdas = self.kin.tilde_lambdas
        
        amp, status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, (0, 1))
        
        if verbose:
            print("\n" + "="*60)
            print("AMPLITUDE COMPUTATION")
            print("="*60)
            
            if amp is not None:
                try:
                    amp_float = float(amp)
                    print(f"  Hodges amplitude: {amp_float:.6e}")
                    print(f"  Sign: {'POSITIVE' if amp_float > 0 else 'NEGATIVE'}")
                except:
                    print(f"  Hodges amplitude: {amp}")
            else:
                print(f"  Failed: {status}")
        
        return amp, status
    
    def canonical_form_ansatz(self, verbose=True):
        """
        Construct an ansatz for the canonical form in kinematic space.
        
        Ω = N / Π_I s_I
        
        where:
        - N is the numerator (may depend on brackets)
        - s_I are the pole invariants
        
        For MHV gravity:
        - Numerator ~ <12>^8 × (something)
        - Denominator = product of all poles
        """
        if verbose:
            print("\n" + "="*60)
            print("CANONICAL FORM ANSATZ")
            print("="*60)
        
        # Helicity factor
        ang_12 = self.kin.angle(0, 1)
        helicity = ang_12**8
        
        # Denominator: product of all poles
        poles = self.factorization_poles()
        
        # Count pole degrees
        # Two-particle poles appear with power 1
        # Three-particle poles also appear with power 1
        
        if verbose:
            print(f"\nNumerator:")
            print(f"  <12>^8 = {float(helicity):.6e}")
            print(f"\nDenominator structure:")
            print(f"  {len([p for p in poles if p['type'] == 'two_particle'])} two-particle poles")
            print(f"  {len([p for p in poles if p['type'] == 'three_particle'])} three-particle poles")
        
        # Compute what the amplitude looks like divided by product of poles
        pole_product = QQ(1)
        for pole in poles:
            if pole['value'] != 0:
                pole_product *= pole['value']
        
        amp, status = self.compute_amplitude(verbose=False)
        
        if amp is not None and pole_product != 0:
            ratio = amp / pole_product
            
            if verbose:
                print(f"\nAmplitude / (product of poles):")
                print(f"  = {float(ratio):.6e}")
                print(f"  (This should relate to the numerator)")
        
        return {
            'helicity_factor': helicity,
            'poles': poles,
            'amplitude': amp
        }
    
    def analyze_positivity(self, verbose=True):
        """
        Analyze positivity structure in kinematic space.
        
        Key question: Is there a region in kinematic space where
        the amplitude is positive (or has definite sign)?
        """
        if verbose:
            print("\n" + "="*60)
            print("POSITIVITY ANALYSIS IN KINEMATIC SPACE")
            print("="*60)
        
        amp, status = self.compute_amplitude(verbose=False)
        boundary = self.boundary_structure(verbose=False)
        
        if amp is None:
            if verbose:
                print("  Cannot analyze: amplitude computation failed")
            return None
        
        try:
            amp_float = float(amp)
        except:
            amp_float = None
        
        # Hypothesis: Amplitude sign correlates with pole signs
        # If all s_ij have same sign as amplitude, we might be in "positive region"
        
        matching_signs = 0
        opposite_signs = 0
        
        amp_sign = 1 if amp_float > 0 else -1 if amp_float < 0 else 0
        
        for pole in boundary['poles']:
            try:
                pole_val = float(pole['value'])
                pole_sign = 1 if pole_val > 0 else -1 if pole_val < 0 else 0
                
                if pole_sign == amp_sign:
                    matching_signs += 1
                elif pole_sign == -amp_sign:
                    opposite_signs += 1
            except:
                pass
        
        if verbose:
            print(f"\nAmplitude sign: {'positive' if amp_sign > 0 else 'negative' if amp_sign < 0 else 'zero'}")
            print(f"Poles with same sign as amplitude: {matching_signs}")
            print(f"Poles with opposite sign: {opposite_signs}")
            
            if matching_signs > opposite_signs:
                print("\n→ Amplitude sign correlates with pole signs!")
            else:
                print("\n→ No clear correlation between amplitude and pole signs")
        
        return {
            'amplitude_sign': amp_sign,
            'matching_signs': matching_signs,
            'opposite_signs': opposite_signs
        }


def run_kinematic_analysis(seeds=[42, 100, 200, 300, 400]):
    """
    Run kinematic space analysis on multiple kinematics.
    """
    print("="*70)
    print("KINEMATIC SPACE GEOMETRY ANALYSIS")
    print("="*70)
    
    results = []
    
    for seed in seeds:
        print(f"\n{'#'*70}")
        print(f"SEED: {seed}")
        print(f"{'#'*70}")
        
        ksg = KinematicSpaceGeometry(seed=seed)
        
        # Run analysis
        ksg.boundary_structure()
        ksg.check_physical_region()
        ksg.compute_amplitude()
        ksg.canonical_form_ansatz()
        positivity = ksg.analyze_positivity()
        
        results.append({
            'seed': seed,
            'positivity': positivity
        })
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    for r in results:
        p = r['positivity']
        if p:
            print(f"  Seed {r['seed']}: amp_sign={p['amplitude_sign']}, matching={p['matching_signs']}, opposite={p['opposite_signs']}")
    
    return results


if __name__ == "__main__":
    run_kinematic_analysis()

