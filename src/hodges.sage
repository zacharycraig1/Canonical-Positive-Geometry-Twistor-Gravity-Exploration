#!/usr/bin/env sage
# =============================================================================
# HODGES MODULE: MomentumTwistor class and Hodges det' formula
# =============================================================================
# Reference: Hodges, arXiv:1204.1930 
# "A simple formula for gravitational MHV amplitudes"
# =============================================================================

from sage.all import *
from itertools import combinations

# Load dependencies
# We use try-except to handle both direct execution and loading
try:
    from src.hodges_sigma import hodges_sigma
except ImportError:
    # Fallback for Sage load() context where src is not a package
    # This assumes CWD is the repo root
    try:
        load('src/hodges_sigma.sage')
    except Exception:
        # If we are in src directory?
        pass

class MomentumTwistor:
    """
    Momentum twistor representation for n particles.
    
    Stores Z_i ∈ QQ^4 for i = 0, ..., n-1.
    Precomputes all angle brackets <i j> and 4-brackets <i j k l>.
    
    Attributes:
        n: Number of particles
        Z: List of momentum twistor vectors
        angle: Dict of angle brackets {(i,j): <i j>}
        four_bracket: Dict of 4-brackets {(i,j,k,l): <i j k l>}
        domain_ok: Whether all required denominators are nonzero
        domain_reason: Explanation if domain_ok is False
    """
    
    def __init__(self, n=6, seed=None, Z=None, check_domain=True):
        """
        Initialize momentum twistor data.
        
        Args:
            n: Number of particles
            seed: Random seed for fallback sampling
            Z: Pre-computed twistor data (list of vectors)
            check_domain: Whether to validate domain
        """
        import numpy as np
        
        self.n = n
        self.check_domain = check_domain
        
        if Z is not None:
            self.Z = Z
        else:
            if seed is not None:
                np.random.seed(seed)
            self.Z = []
            for i in range(n):
                z = vector(QQ, [
                    QQ(np.random.randint(-10, 11)),
                    QQ(np.random.randint(-10, 11)),
                    QQ(np.random.randint(-10, 11)),
                    QQ(np.random.randint(-10, 11))
                ])
                while all(x == 0 for x in z):
                    z = vector(QQ, [
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11))
                    ])
                self.Z.append(z)
        
        self._compute_brackets()
        
        if check_domain:
            self.domain_ok, self.domain_reason = self._check_domain()
        else:
            self.domain_ok = None
            self.domain_reason = None
    
    def _compute_brackets(self):
        """Precompute all brackets once."""
        n = self.n
        
        # Angle brackets: <i j> = Z_i[0]*Z_j[1] - Z_i[1]*Z_j[0]
        self.angle = {}
        for i in range(n):
            for j in range(n):
                self.angle[(i, j)] = self.Z[i][0] * self.Z[j][1] - self.Z[i][1] * self.Z[j][0]
        
        # 4-brackets: <i j k l> = det([Z_i; Z_j; Z_k; Z_l])
        self.four_bracket = {}
        for ijkl in combinations(range(n), 4):
            i, j, k, l = sorted(ijkl)
            M = matrix(QQ, [self.Z[i], self.Z[j], self.Z[k], self.Z[l]])
            self.four_bracket[ijkl] = M.det()
    
    def _check_domain(self):
        """
        Tight domain check: only what formulas actually use.
        
        Checks:
        - All consecutive angle brackets <i, i+1> != 0
        - All square bracket denominators != 0
        """
        n = self.n
        
        # Check consecutive angle brackets
        for i in range(n):
            j = (i + 1) % n
            ang = self.get_angle(i, j)
            if ang == 0:
                return False, f"domain_violation_angle_consecutive_{i}_{j}"
        
        # Check square bracket denominators
        for i in range(n):
            for j in range(n):
                if i != j:
                    im1 = (i - 1) % n
                    jm1 = (j - 1) % n
                    ang_i = self.get_angle(im1, i)
                    ang_j = self.get_angle(jm1, j)
                    if ang_i == 0:
                        return False, f"domain_violation_square_den_{im1}_{i}"
                    if ang_j == 0:
                        return False, f"domain_violation_square_den_{jm1}_{j}"
        
        return True, "ok"
    
    def get_angle(self, i, j):
        """Get angle bracket <i j>."""
        return self.angle.get((i, j), QQ(0))
    
    def get_four_bracket(self, i, j, k, l):
        """
        Get 4-bracket <i j k l> with correct sign.
        
        The sign accounts for the permutation from sorted order.
        """
        indices = tuple(sorted([i, j, k, l]))
        base = self.four_bracket.get(indices, QQ(0))
        
        # Compute sign from permutation
        perm_list = [i, j, k, l]
        sorted_list = sorted(perm_list)
        inversions = 0
        for a in range(len(perm_list)):
            for b in range(a + 1, len(perm_list)):
                if sorted_list.index(perm_list[a]) > sorted_list.index(perm_list[b]):
                    inversions += 1
        sign = 1 if inversions % 2 == 0 else -1
        
        return sign * base
    
    def get_lambda(self, i):
        """Get lambda spinor (top 2 components of Z)."""
        return vector(QQ, [self.Z[i][0], self.Z[i][1]])

    def get_mu(self, i):
        """Get mu spinor (bottom 2 components of Z)."""
        return vector(QQ, [self.Z[i][2], self.Z[i][3]])

    def get_tilde_lambda(self, i):
        """
        Get tilde lambda spinor reconstructed from twistors.
        
        Using formula for standard infinity twistor:
        |i] = ( <i,i+1> mu_{i-1} + <i+1,i-1> mu_i + <i-1,i> mu_{i+1} ) / (<i-1,i><i,i+1>)
        """
        n = self.n
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        
        mu_im1 = self.get_mu(im1)
        mu_i = self.get_mu(i)
        mu_ip1 = self.get_mu(ip1)
        
        ang_i_ip1 = self.get_angle(i, ip1)
        ang_ip1_im1 = self.get_angle(ip1, im1)
        ang_im1_i = self.get_angle(im1, i)
        
        denom = ang_im1_i * ang_i_ip1
        if denom == 0:
            return None
            
        # Numerator vector sum
        num = mu_im1 * ang_i_ip1 + mu_i * ang_ip1_im1 + mu_ip1 * ang_im1_i
        
        return num / denom

    def get_square(self, i, j):
        """
        Get square bracket [i j].
        
        [i j] = det(tilde_lambda_i, tilde_lambda_j)
        """
        lam_t_i = self.get_tilde_lambda(i)
        lam_t_j = self.get_tilde_lambda(j)
        
        if lam_t_i is None or lam_t_j is None:
            return None
            
        # det([u, v]) = u0*v1 - u1*v0
        return lam_t_i[0] * lam_t_j[1] - lam_t_i[1] * lam_t_j[0]


def mandelstam_invariant(twistor, i, j):
    """
    Compute Mandelstam invariant s_{ij}.
    
    s_{ij} = <i j> * [i j]
    
    Args:
        twistor: MomentumTwistor instance
        i, j: Particle indices
        
    Returns:
        s_{ij} as QQ, or None if undefined
    """
    ang_ij = twistor.get_angle(i, j)
    sq_ij = twistor.get_square(i, j)
    if sq_ij is None:
        return None
    return ang_ij * sq_ij


def build_hodges_phi(twistor, ref_spinors):
    """
    Build the full n x n Hodges Phi matrix.
    
    Args:
        twistor: MomentumTwistor instance
        ref_spinors: Tuple (lambda_x, lambda_y) of reference spinors
        
    Returns:
        Phi matrix (n x n) or (None, reason) if error
    """
    n = twistor.n
    Phi = matrix(QQ, n, n)
    lambda_x, lambda_y = ref_spinors
    
    # Extract external spinors λ_i
    lambdas = []
    for i in range(n):
        lambda_i = vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]])
        lambdas.append(lambda_i)
    
    # First, compute all off-diagonal elements
    # Phi_{ij} = [i j] / <i j>
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0:
                    return (None, "domain_violation_angle_bracket_offdiag")
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None:
                    return (None, "domain_violation_square_bracket")
                Phi[i, j] = ij_sq / ij_ang

    
    # Compute diagonal elements using reference spinors λ_x, λ_y (CORRECTED: no special cases)
    # Formula: Φ_{ii} = - Σ_{j≠i} Φ_{ij} * (<j x> <j y>) / (<i x> <i y>)
    for i in range(n):
        lambda_i = lambdas[i]
        ang_ix = ang_vec(lambda_i, lambda_x)  # <i x>
        ang_iy = ang_vec(lambda_i, lambda_y)  # <i y>
        
        if ang_ix == 0 or ang_iy == 0:
            return (None, "domain_violation_angle_bracket_diag")
        
        diag_sum = QQ(0)
        for j in range(n):
            if j == i:
                continue
            lambda_j = lambdas[j]
            ang_jx = ang_vec(lambda_j, lambda_x)  # <j x>
            ang_jy = ang_vec(lambda_j, lambda_y)  # <j y>
            
            if ang_jx == 0 or ang_jy == 0:
                return (None, "domain_violation_angle_bracket_diag")
            
            # Φ_{ii} = - Σ_{j≠i} Φ_{ij} * (<j x> <j y>) / (<i x> <i y>)
            contrib = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= contrib
        
        Phi[i, i] = diag_sum
        
    return Phi

def hodges_6pt_mhv(twistor, deletion_set=None):
    """
    Hodges formula for 6-point MHV gravity - CORRECT det' implementation.
    
    Args:
        twistor: MomentumTwistor instance
        deletion_set: Optional list of 3 indices to delete. Default [0, 1, 2].
        
    Returns:
        Tuple (amplitude, reason) where reason is "ok" or error description
    """
    n = twistor.n
    
    # Reference spinors (2-vectors), not leg indices
    lambda_x, lambda_y = sample_reference_spinors(twistor)
    if lambda_x is None or lambda_y is None:
        return (None, "failed_to_sample_reference_spinors")
        
    Phi = build_hodges_phi(twistor, (lambda_x, lambda_y))
    if isinstance(Phi, tuple): # Error occurred
        return Phi
    
    # Compute reduced determinant det'(Phi)
    if deletion_set is None:
        rows_to_delete = [0, 1, 2]
    else:
        rows_to_delete = deletion_set
        
    cols_to_delete = rows_to_delete # Symmetric deletion for det'
    
    rows_to_keep = [i for i in range(n) if i not in rows_to_delete]
    cols_to_keep = [i for i in range(n) if i not in cols_to_delete]
    
    Phi_red = Phi[rows_to_keep, cols_to_keep]
    
    try:
        det_Phi_red = Phi_red.det()
    except Exception:
        return (None, "determinant_computation_failed")
    
    # Normalization factor: (<ij><jk><ki>)^2 for deleted rows {i,j,k}
    r1, r2, r3 = rows_to_delete
    ang_12 = twistor.get_angle(r1, r2)
    ang_23 = twistor.get_angle(r2, r3)
    ang_31 = twistor.get_angle(r3, r1)
    
    if ang_12 == 0 or ang_23 == 0 or ang_31 == 0:
        return (None, "domain_violation_angle_bracket_normalization")
    norm_factor = (ang_12 * ang_23 * ang_31) ** 2
    
    det_prime_Phi = det_Phi_red / norm_factor
    
    # Final Amplitude Construction
    # M_Grav = det'(Phi) * <0 1>^8  (For MHV with 0,1 negative helicity)
    helicity_factor = twistor.get_angle(0, 1) ** 8
    
    return (det_prime_Phi * helicity_factor, "ok")

def sample_reference_spinors(twistor, max_attempts=100):
    """
    Sample reference spinors λ_x, λ_y such that <i x> and <i y> are nonzero for all i.
    
    Returns: (λ_x, λ_y) as 2-vectors in QQ^2, or (None, None) if failed
    """
    import random
    random.seed(int(42))  # Deterministic for reproducibility
    
    for attempt in range(max_attempts):
        # Sample small integer spinors
        lambda_x = vector(QQ, [QQ(random.randint(1, 10)), QQ(random.randint(1, 10))])
        lambda_y = vector(QQ, [QQ(random.randint(1, 10)), QQ(random.randint(1, 10))])
        
        # Check that all <i x> and <i y> are nonzero
        all_ok = True
        for i in range(twistor.n):
            lambda_i = vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]])
            ang_ix = lambda_i[0] * lambda_x[1] - lambda_i[1] * lambda_x[0]
            ang_iy = lambda_i[0] * lambda_y[1] - lambda_i[1] * lambda_y[0]
            if ang_ix == 0 or ang_iy == 0:
                all_ok = False
                break
        
        if all_ok:
            return (lambda_x, lambda_y)
    
    return (None, None)

def ang_vec(a, b):
    """Angle bracket for 2-vectors: <a b> = a[0]*b[1] - a[1]*b[0]"""
    return a[0] * b[1] - a[1] * b[0]

def hodges_6pt_mhv_reduced(twistor, ref_spinors=None, deletion_set=None):
    """
    Hodges reduced amplitude for 6-point MHV gravity (bar M_6).
    
    Implements Hodges' reduced determinant form:
      \bar M_n = (-1)^{n+1} * σ(ijk,rst) * c_{ijk} * c_{rst} * |Φ|^{rst}_{ijk}
    where:
      c_{ijk} = 1 / (<i j><j k><k i>)
      |Φ|^{rst}_{ijk} is the (n-3)x(n-3) minor after deleting rows i,j,k and columns r,s,t.
    
    Args:
        twistor: MomentumTwistor instance (or adapter)
        ref_spinors: Tuple (lambda_x, lambda_y) of reference spinors. If None, sampled.
        deletion_set: Tuple (rows_delete, cols_delete). If None, uses default ([0,1,2], [3,4,5]).
                     rows_delete: indices of rows to remove (i,j,k)
                     cols_delete: indices of columns to remove (r,s,t)
    
    CORRECTED: Uses reference spinors λ_x, λ_y (2-vectors), not reference leg indices.
    """
    n = twistor.n
    if n != 6:
        return (None, "unsupported_n_for_reduced_hodges")

    # Build full n×n Phi matrix over the same ring as twistor (assumed QQ)
    Phi = matrix(QQ, n, n)

    # Reference spinors (2-vectors), not leg indices
    if ref_spinors is None:
        lambda_x, lambda_y = sample_reference_spinors(twistor)
    else:
        lambda_x, lambda_y = ref_spinors
        
    if lambda_x is None or lambda_y is None:
        return (None, "failed_to_sample_reference_spinors")
    
    # Extract external spinors λ_i
    lambdas = []
    # If twistor has explicit lambdas (Adapter), use them. Else extract from Z.
    if hasattr(twistor, 'lambdas'):
         lambdas = twistor.lambdas
    else:
        for i in range(n):
            lambda_i = vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]])
            lambdas.append(lambda_i)

    # Off-diagonal: Phi_{ij} = [i j] / <i j>
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0:
                    return (None, "domain_violation_angle_bracket_offdiag")
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None:
                    return (None, "domain_violation_square_bracket")
                Phi[i, j] = ij_sq / ij_ang

    # Diagonal using reference spinors λ_x, λ_y (CORRECTED: no special cases)
    # Formula: Φ_{ii} = - Σ_{j≠i} Φ_{ij} * (<j x> <j y>) / (<i x> <i y>)
    for i in range(n):
        lambda_i = lambdas[i]
        ang_ix = ang_vec(lambda_i, lambda_x)  # <i x>
        ang_iy = ang_vec(lambda_i, lambda_y)  # <i y>
        
        if ang_ix == 0 or ang_iy == 0:
            return (None, "domain_violation_angle_bracket_diag")
        
        diag_sum = QQ(0)
        for j in range(n):
            if j == i:
                continue
            lambda_j = lambdas[j]
            ang_jx = ang_vec(lambda_j, lambda_x)  # <j x>
            ang_jy = ang_vec(lambda_j, lambda_y)  # <j y>
            
            if ang_jx == 0 or ang_jy == 0:
                return (None, "domain_violation_angle_bracket_diag")
            
            # Φ_{ii} = - Σ_{j≠i} Φ_{ij} * (<j x> <j y>) / (<i x> <i y>)
            contrib = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= contrib
        
        Phi[i, i] = diag_sum

    # Choose (i,j,k) = rows_delete, (r,s,t) = cols_delete
    if deletion_set is None:
        rows_delete = [0, 1, 2]
        cols_delete = [3, 4, 5]
    else:
        rows_delete, cols_delete = deletion_set
        
    # Determine kept indices
    rows_keep = [r for r in range(n) if r not in rows_delete]
    cols_keep = [c for c in range(n) if c not in cols_delete]
    
    Phi_minor = Phi[rows_keep, cols_keep]
    try:
        det_minor = Phi_minor.det()
    except Exception:
        return (None, "determinant_computation_failed")

    # Compensation factors: c_{ijk} and c_{rst}
    # c_{ijk} = 1 / (<ij><jk><ki>)
    def compute_c_factor(indices):
        i, j, k = indices
        a_ij = twistor.get_angle(i, j)
        a_jk = twistor.get_angle(j, k)
        a_ki = twistor.get_angle(k, i)
        if a_ij == 0 or a_jk == 0 or a_ki == 0:
            return None
        return QQ(1) / (a_ij * a_jk * a_ki)

    cI = compute_c_factor(rows_delete)
    if cI is None: return (None, "domain_violation_angle_cI")
    
    cJ = compute_c_factor(cols_delete)
    if cJ is None: return (None, "domain_violation_angle_cJ")

    # Hodges sigma sign: σ(ijk,rst)
    # Assumes hodges_sigma is available from top-level load/import
    sigma = hodges_sigma(rows_delete, cols_delete, n)

    # Overall sign factor (-1)^{n+1} for Hodges' reduced form
    sign = -1  # since n=6 -> (-1)^{7} = -1

    # Full Hodges reduced formula
    bar_M6 = sign * sigma * cI * cJ * det_minor
    return (bar_M6, "ok")

