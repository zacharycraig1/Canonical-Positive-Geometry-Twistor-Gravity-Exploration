#!/usr/bin/env sage
"""
BCFW Recursion for 6-Point MHV Gravity
======================================

Implements the BCFW recursion relations for MHV gravity amplitudes.

For gravity with shift [1,6⟩ (particles 0 and 5 in 0-indexed):
    M_6 = Σ_{channels} M_L × M_R / P²

where the sum is over factorization channels.

Channels for 6-point:
1. s_01 = 0: M_3(0,1,P) × M_5(-P,2,3,4,5) / s_01
2. s_012 = 0: M_4(0,1,2,P) × M_4(-P,3,4,5) / s_012
3. s_0123 = 0: M_5(0,1,2,3,P) × M_3(-P,4,5) / s_0123

For MHV, the lower-point amplitudes have known forms.

Key reference: Britto-Cachazo-Feng-Witten (BCFW) 2005
"""

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


class BCFWGravity:
    """
    BCFW recursion for MHV gravity amplitudes.
    """
    
    def __init__(self, kinematics, shift=(0, 5)):
        """
        Initialize BCFW calculator.
        
        Args:
            kinematics: SpinorKinematics object
            shift: Tuple (a, b) for BCFW shift [a,b⟩
        """
        self.kin = kinematics
        self.n = 6
        self.shift = shift
        
        # Cache spinor data
        self.lambdas = kinematics.lambdas
        self.tilde_lambdas = kinematics.tilde_lambdas
    
    def angle(self, i, j):
        """Angle bracket ⟨ij⟩."""
        return self.kin.angle(i, j)
    
    def square(self, i, j):
        """Square bracket [ij]."""
        return self.kin.square(i, j)
    
    def s(self, *indices):
        """Mandelstam invariant."""
        if len(indices) == 2:
            return self.kin.s(indices[0], indices[1])
        else:
            # Multi-particle invariant
            total = QQ(0)
            for a, i in enumerate(indices):
                for b in range(a + 1, len(indices)):
                    j = indices[b]
                    total += self.kin.s(i, j)
            return total
    
    def mhv_3(self, particles, neg_hel):
        """
        3-point MHV amplitude.
        
        For 3-point with one negative helicity (a) and two positive (b,c):
            M_3^MHV = ⟨ab⟩^4 ⟨ac⟩^4 / ⟨bc⟩^4
        
        For gravity: the anti-holomorphic 3-point is:
            M_3^{-++(+)} = [bc]^4 / ([ab][ac])  or similar
        
        Actually, for real momenta, 3-point amplitudes vanish.
        In BCFW, we work with complex momenta where they're non-zero.
        
        For MHV gravity:
            M_3^{--+} = ⟨12⟩^8 / (⟨13⟩ ⟨23⟩)^4  (schematically)
        """
        if len(particles) != 3:
            raise ValueError("Need exactly 3 particles")
        
        # For BCFW continuation, 3-point amplitudes are special
        # We use the soft limit formula instead
        return QQ(1)  # Normalized for recursion
    
    def mhv_4(self, particles, neg_hel=(0, 1)):
        """
        4-point MHV gravity amplitude.
        
        M_4^{--++} = ⟨ab⟩^8 / (s × t)
        
        where a,b are negative helicity, s = s_{ab}, t = s_{ac} or s_{bc}.
        
        For particles {p0, p1, p2, p3} with neg_hel specifying which are negative:
        """
        if len(particles) != 4:
            raise ValueError("Need exactly 4 particles")
        
        # Map negative helicity indices to local particles
        neg_local = list(neg_hel)
        pos_local = [i for i in range(4) if i not in neg_local]
        
        a, b = [particles[i] for i in neg_local]
        c, d = [particles[i] for i in pos_local]
        
        # Numerator ⟨ab⟩^8
        ang_ab = self.angle(a, b)
        numer = ang_ab**8
        
        # Denominators: s and t channels
        s = self.s(a, b)  # = s_{cd} by momentum conservation
        t = self.s(a, c)  # = s_{bd}
        
        if s == 0 or t == 0:
            return None
        
        return numer / (s * t)
    
    def mhv_5(self, particles, neg_hel=(0, 1)):
        """
        5-point MHV gravity amplitude.
        
        This is more complex - use Hodges formula on subset.
        For simplicity, use the recursion relation or explicit formula.
        """
        # For a complete implementation, we'd compute this via Hodges
        # For now, return placeholder (this limits accuracy)
        return None  # Not implemented
    
    def bcfw_channels(self):
        """
        Enumerate BCFW factorization channels.
        
        For [0,5⟩ shift, channels are:
        - s_01 = 0 (splits as 3 | 5)
        - s_012 = 0 (splits as 4 | 4)
        - s_0123 = 0 (splits as 5 | 3)
        - s_05 = 0 (splits as 3 | 5 with other grouping)
        - s_015 = 0 (splits as 4 | 4)
        - s_0154 = 0 (splits as 5 | 3)
        
        Actually, for [0,5⟩ shift, we pick channels where momentum flow 
        includes particle 0 but not particle 5.
        """
        n = self.n
        a, b = self.shift
        
        channels = []
        
        # Channel selection: consecutive particles starting from 'a'
        # that don't include 'b'
        
        # 2-particle: s_{a, a+1}
        channels.append({
            'type': '2-particle',
            'left': [a, (a+1) % n],
            'right': [i for i in range(n) if i not in [a, (a+1) % n]],
            'invariant': 's_' + ''.join(str(x) for x in [a, (a+1) % n])
        })
        
        # 3-particle: s_{a, a+1, a+2}
        left_3 = [a, (a+1) % n, (a+2) % n]
        if b not in left_3:
            channels.append({
                'type': '3-particle',
                'left': left_3,
                'right': [i for i in range(n) if i not in left_3],
                'invariant': 's_' + ''.join(str(x) for x in left_3)
            })
        
        # 4-particle: s_{a, a+1, a+2, a+3}
        left_4 = [(a + k) % n for k in range(4)]
        if b not in left_4:
            channels.append({
                'type': '4-particle',
                'left': left_4,
                'right': [i for i in range(n) if i not in left_4],
                'invariant': 's_' + ''.join(str(x) for x in left_4)
            })
        
        # Reverse direction channels (starting from b-1 going down)
        # s_{b-1, b} channel
        bm1 = (b - 1) % n
        channels.append({
            'type': '2-particle-rev',
            'left': [bm1, b],
            'right': [i for i in range(n) if i not in [bm1, b]],
            'invariant': 's_' + ''.join(str(x) for x in sorted([bm1, b]))
        })
        
        return channels
    
    def channel_contribution(self, channel, verbose=False):
        """
        Compute the contribution from a single BCFW channel.
        
        Each channel contributes:
            M_left × M_right / P²
        
        where P² is the propagator momentum squared (= invariant for that channel).
        """
        left = channel['left']
        right = channel['right']
        n_left = len(left)
        n_right = len(right)
        
        # Propagator
        P_sq = self.s(*left)
        if P_sq == 0:
            if verbose:
                print(f"  {channel['invariant']}: P² = 0, skip")
            return None
        
        # Compute left amplitude
        if n_left == 3:
            M_left = self.mhv_3(left, neg_hel=(0,))  # Placeholder
        elif n_left == 4:
            # Determine which particles have negative helicity in the subset
            # Original neg hel: particles 0, 1
            neg_in_left = [p for p in left if p in [0, 1]]
            if len(neg_in_left) < 2:
                # Not MHV on left - complex structure
                M_left = QQ(1)  # Placeholder
            else:
                M_left = self.mhv_4(left, neg_hel=tuple(range(len(neg_in_left))))
        elif n_left == 5:
            M_left = self.mhv_5(left)
        else:
            M_left = QQ(1)
        
        # Compute right amplitude (similar logic)
        if n_right == 3:
            M_right = self.mhv_3(right, neg_hel=())
        elif n_right == 4:
            neg_in_right = [p for p in right if p in [0, 1]]
            if len(neg_in_right) < 2:
                M_right = QQ(1)
            else:
                M_right = self.mhv_4(right, neg_hel=tuple(range(len(neg_in_right))))
        elif n_right == 5:
            M_right = self.mhv_5(right)
        else:
            M_right = QQ(1)
        
        if M_left is None or M_right is None:
            if verbose:
                print(f"  {channel['invariant']}: amplitude computation failed")
            return None
        
        contrib = M_left * M_right / P_sq
        
        if verbose:
            try:
                print(f"  {channel['invariant']}: M_L={float(M_left):.4e}, M_R={float(M_right):.4e}, P²={float(P_sq):.4e}")
                print(f"    Contribution: {float(contrib):.6e}")
            except:
                print(f"  {channel['invariant']}: contribution = {contrib}")
        
        return contrib
    
    def compute_bcfw_amplitude(self, verbose=True):
        """
        Compute the 6-point MHV gravity amplitude via BCFW.
        
        M_6 = Σ_{channels} M_L × M_R / P²
        """
        if verbose:
            print("="*60)
            print("BCFW GRAVITY AMPLITUDE")
            print("="*60)
            print(f"Shift: [{self.shift[0]}, {self.shift[1]}⟩")
        
        channels = self.bcfw_channels()
        
        if verbose:
            print(f"\n{len(channels)} BCFW channels:")
            for ch in channels:
                print(f"  {ch['invariant']}: left={ch['left']}, right={ch['right']}")
        
        total = QQ(0)
        contributions = []
        
        if verbose:
            print(f"\nChannel contributions:")
        
        for ch in channels:
            contrib = self.channel_contribution(ch, verbose=verbose)
            if contrib is not None:
                total += contrib
                contributions.append({
                    'channel': ch,
                    'contribution': contrib
                })
        
        # Include helicity factor for proper normalization
        ang_01 = self.angle(0, 1)
        helicity = ang_01**8
        
        if verbose:
            print(f"\nRaw BCFW sum: {float(total):.6e}")
            print(f"Helicity factor ⟨01⟩^8: {float(helicity):.6e}")
        
        # The BCFW result should be compared with Hodges
        return {
            'raw_sum': total,
            'helicity_factor': helicity,
            'channels': contributions
        }
    
    def compare_with_hodges(self, verbose=True):
        """
        Compare BCFW result with Hodges formula.
        """
        if verbose:
            print("\n" + "="*60)
            print("COMPARISON WITH HODGES")
            print("="*60)
        
        # Compute Hodges
        hodges_amp, hodges_status = hodges_npt_mhv_canonical(
            self.lambdas, self.tilde_lambdas, (0, 1)
        )
        
        if verbose:
            if hodges_amp is not None:
                print(f"Hodges amplitude: {float(hodges_amp):.6e}")
            else:
                print(f"Hodges failed: {hodges_status}")
        
        # Compute BCFW
        bcfw_result = self.compute_bcfw_amplitude(verbose=False)
        bcfw_raw = bcfw_result['raw_sum']
        
        if verbose:
            print(f"BCFW raw sum:    {float(bcfw_raw):.6e}")
        
        if hodges_amp is not None and bcfw_raw is not None:
            ratio = bcfw_raw / hodges_amp if hodges_amp != 0 else float('inf')
            if verbose:
                print(f"Ratio (BCFW/Hodges): {float(ratio):.6f}")
        
        return {
            'hodges': hodges_amp,
            'bcfw': bcfw_raw,
            'match': False  # Need to implement proper BCFW formula
        }


def run_bcfw_test(seeds=[42, 100, 200]):
    """
    Test BCFW recursion on multiple kinematics.
    """
    print("="*70)
    print("BCFW GRAVITY TEST")
    print("="*70)
    
    for seed in seeds:
        print(f"\n{'#'*70}")
        print(f"SEED: {seed}")
        print(f"{'#'*70}")
        
        kin = SpinorKinematics.random_rational(6, seed=seed)
        bcfw = BCFWGravity(kin)
        
        result = bcfw.compare_with_hodges()
    
    return True


if __name__ == "__main__":
    run_bcfw_test()

