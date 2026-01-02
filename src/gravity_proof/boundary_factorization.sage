#!/usr/bin/env sage
# =============================================================================
# PHASE 6: Boundary Factorization Verification
# =============================================================================
# Verifies that residues of Omega(R6) at boundaries factorize correctly

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

# Load Sage modules
load("src/gravity_proof/canonical_form_gravity.sage")

class BoundaryFactorization:
    """
    Verifies that boundary residues of canonical form factorize correctly.
    
    At a factorization channel s_{ijk} → 0, the amplitude factorizes:
    M_6 → M_4(i,j,k,P) × 1/s_{ijk} × M_4(-P, remaining)
    
    The geometric interpretation:
    Res_{s_{ijk}=0} Ω(R) = Ω(R_L) × Ω(R_R)
    """
    
    def __init__(self, canonical_form):
        """
        Initialize factorization checker.
        
        Args:
            canonical_form: CanonicalFormR6 instance
        """
        self.canon_form = canonical_form
        self.region = canonical_form.region
        self.kin = canonical_form.kin
    
    def four_point_gravity_mhv(self, particles):
        """
        Compute 4-point MHV gravity amplitude.
        
        Formula from directive line 564:
        M_4^{MHV} = ⟨12⟩^4 [34]^4 / (s_{12} s_{23})
        
        Args:
            particles: tuple (i1, i2, i3, i4) of particle indices (0-indexed)
                      where i1, i2 are negative helicity
        
        Returns:
            M_4^MHV value
        """
        if len(particles) != 4:
            raise ValueError("4-point amplitude requires 4 particles")
        
        i1, i2, i3, i4 = particles
        
        # Negative helicity: i1, i2
        # Positive helicity: i3, i4
        
        # Compute brackets
        ang_12 = self.kin.angle(i1, i2)
        sq_34 = self.kin.square(i3, i4)
        
        # Compute Mandelstam invariants
        s_12 = self.kin.s(i1, i2)
        s_23 = self.kin.s(i2, i3)
        
        # M_4 = ⟨12⟩^4 [34]^4 / (s_{12} s_{23})
        numerator = (ang_12**4) * (sq_34**4)
        denominator = s_12 * s_23
        
        if denominator == 0:
            return None  # On pole
        
        return numerator / denominator
    
    def factorization_channels(self):
        """
        Return list of factorization channels for 6-point.
        
        Channels: s_{123}, s_{234}, s_{345}, s_{456}, s_{156}, s_{126}
        (and cyclic permutations)
        
        Returns:
            List of dicts with channel information
        """
        channels = []
        
        # Standard channels
        channel_specs = [
            {'name': 's123', 'particles': (0, 1, 2), 'remaining': (3, 4, 5)},
            {'name': 's234', 'particles': (1, 2, 3), 'remaining': (0, 4, 5)},
            {'name': 's345', 'particles': (2, 3, 4), 'remaining': (0, 1, 5)},
            {'name': 's456', 'particles': (3, 4, 5), 'remaining': (0, 1, 2)},
            {'name': 's156', 'particles': (0, 4, 5), 'remaining': (1, 2, 3)},
            {'name': 's126', 'particles': (0, 1, 5), 'remaining': (2, 3, 4)},
        ]
        
        for spec in channel_specs:
            # Compute s value
            i, j, k = spec['particles']
            s_val = self.kin.s(i, j) + self.kin.s(i, k) + self.kin.s(j, k)
            
            channels.append({
                'name': spec['name'],
                'particles': spec['particles'],
                'remaining': spec['remaining'],
                's_value': s_val,
                'description': f'Channel where {spec["name"]} → 0'
            })
        
        return channels
    
    def compute_residue_at_channel(self, channel_name):
        """
        Compute residue of canonical form at a factorization channel.
        
        Args:
            channel_name: name of channel (e.g., 's123')
            
        Returns:
            Residue value (should equal product of 4-point amplitudes)
        """
        channels = self.factorization_channels()
        channel = None
        
        for c in channels:
            if c['name'] == channel_name:
                channel = c
                break
        
        if channel is None:
            raise ValueError(f"Channel {channel_name} not found")
        
        # The residue is computed by taking limit as s_{ijk} → 0
        # This is complex - requires understanding how s_{ijk} relates to z coordinates
        
        # For now, return structure
        return {
            'channel': channel_name,
            'residue_value': 'computed_via_limit',
            'note': 'Residue computation requires limit-taking as s → 0'
        }
    
    def verify_factorization(self, channel_name):
        """
        Verify that residue at channel equals product of 4-point amplitudes.
        
        Args:
            channel_name: name of channel
            
        Returns:
            dict with verification results
        """
        channel = None
        for c in self.factorization_channels():
            if c['name'] == channel_name:
                channel = c
                break
        
        if channel is None:
            raise ValueError(f"Channel {channel_name} not found")
        
        # Compute residue
        residue = self.compute_residue_at_channel(channel_name)
        
        # Compute 4-point amplitudes
        # For channel s_{ijk}, we need M_4(i,j,k,P) and M_4(-P, remaining)
        # But P is the intermediate momentum, so we need to construct
        # the 4-point kinematics from the 6-point kinematics
        
        # Left amplitude: particles in channel
        left_particles = channel['particles']
        # Right amplitude: remaining particles
        
        # This is complex - requires mapping from 6-point to 4-point kinematics
        # For now, return structure
        
        return {
            'channel': channel_name,
            'residue': residue,
            'left_amplitude': 'M_4(channel_particles)',
            'right_amplitude': 'M_4(remaining_particles)',
            'factorization_verified': False,
            'note': 'Requires mapping from 6-point to 4-point kinematics at factorization limit'
        }
    
    def verify_all_channels(self):
        """
        Verify factorization for all channels.
        
        Returns:
            dict with results for each channel
        """
        channels = self.factorization_channels()
        results = {}
        
        for channel in channels:
            name = channel['name']
            try:
                result = self.verify_factorization(name)
                results[name] = result
            except Exception as e:
                results[name] = {
                    'channel': name,
                    'error': str(e)
                }
        
        return results


# Test function
def test_boundary_factorization():
    """Quick test of BoundaryFactorization"""
    print("Testing BoundaryFactorization...")
    
    from src.kinematics.spinors import SpinorKinematics
    from src.gravity_proof.psi_matrix import PsiMatrixMHV
    from src.gravity_proof.positivity_region import PositivityRegionR6
    from src.gravity_proof.canonical_form_gravity import CanonicalFormR6
    
    # Create test setup
    kin = SpinorKinematics.random_rational(6, seed=42)
    var('z4 z5 z6')
    psi = PsiMatrixMHV(kin, z4, z5, z6)
    region = PositivityRegionR6(psi)
    canon = CanonicalFormR6(region)
    
    # Create factorization checker
    factor = BoundaryFactorization(canon)
    
    # Test 4-point amplitude
    try:
        m4 = factor.four_point_gravity_mhv((0, 1, 2, 3))
        print(f"4-point MHV amplitude: {m4}")
    except Exception as e:
        print(f"4-point amplitude computation failed: {e}")
    
    # Test channels
    channels = factor.factorization_channels()
    print(f"Number of factorization channels: {len(channels)}")
    for c in channels[:3]:  # Print first 3
        print(f"Channel {c['name']}: {c['description']}")
    
    # Test verification (will likely fail as it's not fully implemented)
    try:
        results = factor.verify_all_channels()
        print(f"Verification results for {len(results)} channels")
    except Exception as e:
        print(f"Verification failed: {e}")
    
    print("Test complete.")


# Only run test if executed directly (commented out to avoid running on load)
# if __name__ == "__main__":
#     test_boundary_factorization()

