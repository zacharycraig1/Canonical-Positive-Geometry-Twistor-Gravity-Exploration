from sage.all import *
from itertools import combinations

def ang_bracket(lambdas, i, j):
    """Compute <i j> = det(lambda_i, lambda_j)."""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]

def sq_bracket(tilde_lambdas, i, j):
    """
    Compute [i j] = det(tilde_lambda_i, tilde_lambda_j).
    Convention: Must match s_ij = <i j>[j i] or <i j>[i j]?
    Standard: s_ij = (p_i + p_j)^2 = 2 p_i.p_j = <i j>[j i].
    Wait, usually s_ij = <ij>[ji].
    But some conventions use s_ij = <ij>[ij] if [ij] is defined as -[ji].
    Let's check src/hodges.sage mandelstam_invariant:
      s_{ij} = <i j> * [i j]
    This implies [i j] in hodges.sage includes the sign flip or convention match.
    
    Standard spinor helicity:
    p_i^{\dot{\alpha}\alpha} = \tilde{\lambda}_i^{\dot{\alpha}} \lambda_i^{\alpha}
    2 p_i \cdot p_j = \langle i j \rangle [j i]
    
    If hodges.sage uses s = <ij>[ij], then its [ij] must be [ji]_{std}.
    Let's stick to calculating [i j]_{std} and then adapting to what hodges expects.
    
    Actually, let's implement [i j] such that s_ij = <i j>[i j] is the standard Mandelstam.
    If s_ij = <ij>[ji], then we define our [ij] function to return [ji]_{std}.
    
    Let's use the explicit determinant for [i j]:
    [i j] = tilde_lambda_i[0]*tilde_lambda_j[1] - tilde_lambda_i[1]*tilde_lambda_j[0]
          = det(tilde_i, tilde_j)
    
    Then <i j>[i j] = det(lam_i, lam_j) * det(til_i, til_j).
    This is antisymmetric in i,j.
    s_ij must be symmetric.
    So s_ij = <i j>[i j] cannot be correct if both are antisymmetric determinants.
    (-1) * (-1) = 1. So <ji>[ji] = <ij>[ij].
    Yes, s_ij is symmetric.
    
    So s_ij = <i j> [i j] is consistent with symmetry.
    Let's assume [i j] = det(tilde_i, tilde_j).
    """
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]

def spinors_to_sij(lambdas, tilde_lambdas):
    """
    Convert spinors to Mandelstam invariants s_ij.
    s_ij = <i j> [i j] (using determinant definition for both).
    """
    n = len(lambdas)
    s_ij = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                s_ij[(i, j)] = QQ(0)
            else:
                ang = ang_bracket(lambdas, i, j)
                sq = sq_bracket(tilde_lambdas, i, j)
                s_ij[(i, j)] = ang * sq
    return s_ij

def sij_to_channels(sij, n=6):
    """
    Compute channel invariants s_S for |S|=2,3.
    s_S = (sum_{i in S} p_i)^2 = sum_{i<j in S} s_ij
    """
    channels = {}
    
    # 2-particle channels (same as s_ij)
    for i in range(n):
        for j in range(i + 1, n):
            # Key is tuple sorted
            key = tuple(sorted((i, j)))
            channels[key] = sij[(i, j)]
            
    # 3-particle channels
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                key = tuple(sorted((i, j, k)))
                # s_ijk = s_ij + s_ik + s_jk
                val = sij[(i, j)] + sij[(i, k)] + sij[(j, k)]
                channels[key] = val
                
    return channels

def spinors_to_channels(lambdas, tilde_lambdas):
    sij = spinors_to_sij(lambdas, tilde_lambdas)
    return sij_to_channels(sij, len(lambdas))

class SpinorHelicityAdapter:
    """
    Adapts lambdas/tilde_lambdas to the interface expected by Hodges/KLT.
    Mimics MomentumTwistor.
    """
    def __init__(self, lambdas, tilde_lambdas):
        self.n = len(lambdas)
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas
        self.Z = [] # Mock Z to allow index access if needed, though we override getters
        # Create mock Z structure to satisfy some checks if they access Z directly
        for i in range(self.n):
            # Z[i] = (lambda_i, mu_i). We only populate lambda part correctly.
            self.Z.append(vector(QQ, [lambdas[i][0], lambdas[i][1], 0, 0]))

    def get_angle(self, i, j):
        return ang_bracket(self.lambdas, i, j)

    def get_square(self, i, j):
        return sq_bracket(self.tilde_lambdas, i, j)
        
    def get_four_bracket(self, i, j, k, l):
        # Not needed for reduced Hodges if we use the one based on [ij]/<ij>
        # But MomentumTwistor computes [ij] from four brackets.
        # Since we override get_square, this might not be called.
        raise NotImplementedError("Four brackets not implemented for direct spinor input")










