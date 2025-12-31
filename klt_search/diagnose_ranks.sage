
from sage.all import *
import numpy as np
import os
import sys

# Add root to path
sys.path.append(os.getcwd())

# Load dependencies
load("src/hodges.sage")
from itertools import permutations

def ts():
    import time
    return time.strftime("%H:%M:%S")

def log_message(msg):
    print(f"[{ts()}] {msg}")

# =============================================================================
# KLT KERNEL (Copied from correct_klt_proof.sage for isolation)
# =============================================================================

def klt_momentum_kernel_6pt(alpha, beta, twistor):
    """
    Compute KLT momentum kernel S_KLT[alpha|beta] for 6-point.
    
    For n=6:
    - Permuted set is {2,3,4} (0-based: {1,2,3})
    - Fixed legs are {1,5,6} (0-based: {0,4,5})
    - alpha, beta are permutations of {1,2,3}
    
    Formula:
    S_KLT[alpha|beta] = 
        (s_{0,alpha[0]} + theta*s_{alpha[0],alpha[1]} + theta*s_{alpha[0],alpha[2]})
      * (s_{0,alpha[1]} + theta*s_{alpha[1],alpha[2]})
      * (s_{0,alpha[2]})
    
    where theta_beta(a,b) = 1 if a appears after b in beta, else 0
    """
    if len(alpha) != 3 or len(beta) != 3:
        return None
    
    # Build position map for beta
    pos_in_beta = {}
    for idx, val in enumerate(beta):
        pos_in_beta[val] = idx
    
    # Theta function: theta_beta(a,b) = 1 if a appears after b in beta
    def theta_beta(a, b):
        if a not in pos_in_beta or b not in pos_in_beta:
            return 0
        return 1 if pos_in_beta[a] > pos_in_beta[b] else 0
    
    # Compute the three factors
    # Factor 1: s_{0,alpha[0]} + theta*s_{alpha[0],alpha[1]} + theta*s_{alpha[0],alpha[2]}
    s_0_a0 = mandelstam_invariant(twistor, 0, alpha[0])
    if s_0_a0 is None: return None
    
    sum1 = s_0_a0
    
    if theta_beta(alpha[0], alpha[1]):
        s_a0_a1 = mandelstam_invariant(twistor, alpha[0], alpha[1])
        if s_a0_a1 is None: return None
        sum1 += s_a0_a1
    
    if theta_beta(alpha[0], alpha[2]):
        s_a0_a2 = mandelstam_invariant(twistor, alpha[0], alpha[2])
        if s_a0_a2 is None: return None
        sum1 += s_a0_a2
    
    # Factor 2: s_{0,alpha[1]} + theta*s_{alpha[1],alpha[2]}
    s_0_a1 = mandelstam_invariant(twistor, 0, alpha[1])
    if s_0_a1 is None: return None
    
    sum2 = s_0_a1
    
    if theta_beta(alpha[1], alpha[2]):
        s_a1_a2 = mandelstam_invariant(twistor, alpha[1], alpha[2])
        if s_a1_a2 is None: return None
        sum2 += s_a1_a2
    
    # Factor 3: s_{0,alpha[2]}
    s_0_a2 = mandelstam_invariant(twistor, 0, alpha[2])
    if s_0_a2 is None: return None
    
    # Product
    return sum1 * sum2 * s_0_a2

# =============================================================================
# GENERICITY CHECKS
# =============================================================================

def check_genericity(twistor):
    """
    Check if the kinematic point is generic.
    1. All s_ij != 0
    2. Gram determinant of 4 independent momenta != 0
    """
    n = twistor.n
    
    # 1. Check all s_ij
    s_values = []
    for i in range(n):
        for j in range(i+1, n):
            s = mandelstam_invariant(twistor, i, j)
            if s is None or s == 0:
                return False, f"s_{i}_{j} is zero or undefined"
            s_values.append(s)
            
    # Check for too many repeated values (heuristic for low-dim kinematics)
    if len(set(s_values)) < 5:
        return False, "Too few unique s_ij values"
        
    # 2. Check Gram determinant
    # Pick p0, p1, p2, p3
    # Gram matrix G_ij = p_i . p_j = s_ij / 2 (assuming massless)
    # Actually p_i . p_j = s_ij / 2
    # But s_ii = 0
    
    # We'll build the Gram matrix for indices 0,1,2,3
    indices = [0, 1, 2, 3]
    G = matrix(QQ, 4, 4)
    for r, i in enumerate(indices):
        for c, j in enumerate(indices):
            if i == j:
                G[r, c] = 0
            else:
                s = mandelstam_invariant(twistor, i, j)
                if s is None: return False, "undefined s_ij in Gram"
                G[r, c] = s / 2
                
    gdet = G.det()
    if gdet == 0:
        return False, "Gram determinant is 0"
        
    return True, "ok"

# =============================================================================
# DIRECT SPINOR SAMPLER (Bypass Twistors)
# =============================================================================

class DirectSpinors:
    def __init__(self, n=6, seed=None):
        self.n = n
        if seed is not None:
            np.random.seed(seed)
            set_random_seed(seed)
            
        # 1. Generate random Lambda (2 x n)
        # We use integers for exactness if possible, but null space might introduce rationals
        self.Lam = matrix(QQ, 2, n, [np.random.randint(-10, 11) for _ in range(2*n)])
        
        # Ensure Lambda is rank 2
        while self.Lam.rank() < 2:
            self.Lam = matrix(QQ, 2, n, [np.random.randint(-10, 11) for _ in range(2*n)])
            
        # 2. Find Null Space (Ker(Lam)) -> dimension n-2 = 4
        # right_kernel() returns vectors v such that Lam * v = 0
        K = self.Lam.right_kernel()
        basis = K.basis() # List of n-vectors
        
        # 3. Generate random TildeLambda (2 x n) from Null Space
        # We need 2 rows, each in the null space of Lambda.
        # Actually, momentum conservation: sum |i> [i| = 0
        # This implies Lam * LamTilde^T = 0
        # So rows of LamTilde must be in the kernel of Lam (viewed appropriately)? 
        # Wait. p_i = |i>[i|. sum p_i = 0.
        # sum lam_i alpha * lam_tilde_i dot_alpha = 0
        # Matrix form: Lam * LamTilde^T = 0 (2x2 zero matrix)
        # Yes. So columns of LamTilde^T (which are rows of LamTilde) must be in the kernel of Lam.
        
        # We construct LamTilde^T (n x 2) by taking linear combos of basis vectors
        # coeffs = matrix(QQ, n-2, 2)
        coeffs = matrix(QQ, len(basis), 2, [np.random.randint(-10, 11) for _ in range(len(basis)*2)])
        
        # Basis matrix B is (n-2) x n
        B = matrix(QQ, basis) 
        
        # LamTilde^T = B.T * coeffs ? No.
        # v in Kernel means Lam * v = 0.
        # We want cols of LamTilde^T to be v.
        # So LamTilde^T = Sum c_k * v_k
        # Let's verify dimensions.
        # v_k is length n.
        # LamTilde^T is n x 2.
        # Col 0 = sum a_k v_k
        # Col 1 = sum b_k v_k
        
        LT_T = B.transpose() * coeffs # (n x n-2) * (n-2 x 2) = n x 2
        self.LamTilde = LT_T.transpose() # 2 x n
        
        # Precompute brackets
        self._compute_brackets()
        self.Z = [] # Mock Z to allow accessing components if needed, but we rely on brackets
        # Fill Z with mock data to satisfy checks if they inspect Z directly
        for i in range(n):
            self.Z.append(vector(QQ, [self.Lam[0,i], self.Lam[1,i], 0, 0]))

    def _compute_brackets(self):
        self.angle = {}
        self.square = {}
        
        for i in range(self.n):
            for j in range(self.n):
                # <i j> = det(lam_i, lam_j)
                val_ang = self.Lam[0, i] * self.Lam[1, j] - self.Lam[1, i] * self.Lam[0, j]
                self.angle[(i, j)] = val_ang
                
                # [i j] = det(lam_tilde_i, lam_tilde_j)
                # Note: Convention signs matter. [i j] = tlam_i_1 * tlam_j_2 - ...
                # Usually defined as det.
                val_sq = self.LamTilde[0, i] * self.LamTilde[1, j] - self.LamTilde[1, i] * self.LamTilde[0, j]
                self.square[(i, j)] = val_sq

    def get_angle(self, i, j):
        return self.angle.get((i, j), 0)
        
    def get_square(self, i, j):
        return self.square.get((i, j), 0)
    
    # Mocking MomentumTwistor interface
    @property
    def lambdas(self):
        # Return list of vectors for reference spinor sampling
        lams = []
        for i in range(self.n):
            lams.append(vector(QQ, [self.Lam[0,i], self.Lam[1,i]]))
        return lams

# Override mandelstam_invariant to use get_square directly without checking for None
# But the imported one might check for None. 
# The imported mandelstam_invariant uses twistor.get_square.
# Our get_square returns 0 if missing, which is fine (or should it be None?)
# The existing one returns None if denominator is zero. Here we have no denominator.

# =============================================================================
# DIAGNOSIS
# =============================================================================

def diagnose_ranks(n_samples=20):
    log_message(f"Starting diagnosis with {n_samples} generic samples...")
    
    passed_phi = 0
    passed_klt = 0
    passed_det_indep = 0
    
    for i in range(n_samples):
        # USE FIXED MOMENTUM TWISTOR
        twistor = MomentumTwistor(n=6, seed=i+2000)
        
        # Check genericity
        is_generic, reason = check_genericity(twistor)
        if not is_generic:
            log_message(f"Sample {i}: SKIPPED - {reason}")
            continue
            
        # 1. Check Hodges Phi Rank
        ref_spinors = sample_reference_spinors(twistor)
        if ref_spinors[0] is None:
            log_message(f"Sample {i}: Failed to sample reference spinors")
            continue
            
        Phi = build_hodges_phi(twistor, ref_spinors)
        if isinstance(Phi, tuple):
            log_message(f"Sample {i}: Failed to build Phi - {Phi[1]}")
            continue
            
        rk_phi = Phi.rank()
        if rk_phi == 3:
            passed_phi += 1
        else:
            log_message(f"Sample {i}: Phi rank = {rk_phi} (EXPECTED 3)")
            
        # 2. Check Deletion Independence
        # det' = det(minor) / norm_factor
        # Check two different deletions: (0,1,2) and (3,4,5)
        
        # Deletion 1: (0,1,2)
        rows1, cols1 = [3,4,5], [3,4,5]
        det1 = Phi[rows1, cols1].det()
        # Norm factor for (0,1,2): (<01><12><20>)^2
        n1 = (twistor.get_angle(0,1) * twistor.get_angle(1,2) * twistor.get_angle(2,0))**2
        val1 = det1 / n1
        
        # Deletion 2: (1,2,3) -> rows [0,4,5]
        rows2, cols2 = [0,4,5], [0,4,5] # Keeping 0,4,5
        det2 = Phi[rows2, cols2].det()
        # Norm factor for (1,2,3): (<12><23><31>)^2
        n2 = (twistor.get_angle(1,2) * twistor.get_angle(2,3) * twistor.get_angle(3,1))**2
        val2 = det2 / n2
        
        # Check equality
        ratio_det = val1 / val2
        if ratio_det == 1:
            passed_det_indep += 1
        else:
             log_message(f"Sample {i}: Det' independence FAILED. Ratio = {float(ratio_det)}")

        # 3. Check KLT Kernel Rank
        from itertools import permutations
        permuted_set = [1, 2, 3] # indices
        perms = sorted(list(permutations(permuted_set)))
        dim_klt = len(perms)
        S_mat = matrix(QQ, dim_klt, dim_klt)
        
        for r, alpha in enumerate(perms):
            for c, beta in enumerate(perms):
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), twistor)
                S_mat[r, c] = val
                
        rk_klt = S_mat.rank()
        if rk_klt == 6:
            passed_klt += 1
        else:
            log_message(f"Sample {i}: KLT Kernel rank = {rk_klt} (EXPECTED 6)")

    log_message("-" * 40)
    log_message(f"Phi Rank 3: {passed_phi}/{n_samples}")
    log_message(f"Det' Independent: {passed_det_indep}/{n_samples}")
    log_message(f"KLT Rank 6: {passed_klt}/{n_samples}")

def verify_amplitude_ratio(n_samples=20):
    log_message(f"\nVerifying KLT/Hodges Amplitude Ratio on {n_samples} points...")
    
    # We need the full KLT amplitude function. 
    # Importing it here or defining it. 
    # For now, let's pull it from correct_klt_proof.sage if possible, or redefine simplified version.
    
    # We'll use the definition:
    # M_KLT = sum_{alpha, beta} PT(1,alpha,n,n-1) * S[alpha|beta] * PT(1,beta,n-1,n)
    # Note: KLT usually uses PT(..., n-1, n) and PT(..., n, n-1) or similar.
    # The one in correct_klt_proof.sage was:
    # M_6 = Σ A(1,alpha,5,6) * S * A(1,beta,6,5)
    # where alpha,beta on {2,3,4}. Fixed {1,5,6}.
    
    # from src.hodges import hodges_6pt_mhv  <-- REMOVED
    from itertools import permutations
    
    ratios = []
    
    for i in range(n_samples):
        twistor = MomentumTwistor(n=6, seed=i+3000)
        is_generic, _ = check_genericity(twistor)
        if not is_generic: continue
        
        # 1. Hodges
        # Use deletion [0,1,2]
        H_res = hodges_6pt_mhv(twistor)
        if isinstance(H_res, tuple): H = H_res[0]
        else: H = H_res
        
        if H is None or H == 0: continue
        
        # 2. KLT
        # M_6 = sum_{alpha,beta} PT(0, alpha, 4, 5) * S[alpha|beta] * PT(0, beta, 5, 4)
        # alpha, beta perm of {1,2,3}
        
        # PT function
        def get_pt(order):
            # A_MHV = <0 1>^4 / prod <i i+1>
            # Assuming 0,1 are neg helicity.
            # Denom is cyclic product of order.
            denom = QQ(1)
            for k in range(6):
                denom *= twistor.get_angle(order[k], order[(k+1)%6])
            if denom == 0: return None
            # MHV numerator: <0 1>^4
            num = twistor.get_angle(0, 1)**4
            return num / denom

        permuted = [1, 2, 3]
        perms = list(permutations(permuted))
        
        M_KLT = QQ(0)
        for alpha in perms:
            # A_L = PT(0, alpha, 4, 5)
            order_L = [0] + list(alpha) + [4, 5]
            AL = get_pt(order_L)
            if AL is None: continue
            
            for beta in perms:
                # A_R = PT(0, beta, 5, 4)  <-- NOTE SWAP 5,4
                order_R = [0] + list(beta) + [5, 4]
                AR = get_pt(order_R)
                if AR is None: continue
                
                S = klt_momentum_kernel_6pt(list(alpha), list(beta), twistor)
                if S is None: continue
                
                M_KLT += AL * S * AR
                
        if M_KLT == 0: continue
        
        ratio = M_KLT / H
        ratios.append(ratio)
        log_message(f"Sample {i}: Ratio = {float(ratio):.6f} (Exact: {ratio})")
        
    if len(ratios) > 0:
        unique_ratios = set(ratios)
        if len(unique_ratios) == 1:
            log_message(f"\n[SUCCESS] Constant Ratio Found: {ratios[0]}")
        else:
             log_message(f"\n[WARNING] Ratios vary. Unique values: {len(unique_ratios)}")
             # check variance
             first = ratios[0]
             if all(abs((r - first)/first) < 1e-10 for r in ratios):
                 log_message("  (But numerically close - likely precision/normalization issue if not exact)")

if __name__ == "__main__":
    diagnose_ranks()
    # verify_amplitude_ratio()
    
    # Test Moment Curve Points (Cleaner)
    log_message("\nTesting Moment Curve Points (Positive Geometry Candidates)...")
    # from src.hodges import sample_reference_spinors  <-- REMOVED
    
    # Moment curve sampler locally to avoid dependency issues
    def get_moment_curve(seed):
        # t = [1, ..., 6] + perturbation
        import numpy as np
        np.random.seed(seed)
        ts = [i + 1 + np.random.normal(0, 0.01) for i in range(6)]
        ts.sort()
        Z = []
        for t in ts:
            Z.append(vector(QQ, [1, t, t**2, t**3])) # Actually this is P^3, we need twistors in P^3?
            # Standard moment curve in P3 is (1, t, t^2, t^3).
            # This gives 4-vectors. Matches momentum twistor format.
        return MomentumTwistor(n=6, Z=Z, check_domain=False) # Skip strict domain check for floats? No, convert to QQ.
        
    def get_moment_curve_rational(seed):
        t = [QQ(i+1) for i in range(6)]
        # Perturb
        import random
        random.seed(seed)
        t = [x + QQ(random.randint(1,100))/1000 for x in t]
        t.sort()
        Z = []
        for val in t:
             Z.append(vector(QQ, [1, val, val**2, val**3]))
        return MomentumTwistor(n=6, Z=Z)

    # Use existing verify_amplitude_ratio logic but with moment curve
    log_message("Using Rational Moment Curve...")
    
    ratios_mc = []
    for i in range(10):
        twistor = get_moment_curve_rational(i)
        
        # Hodges
        H_res = hodges_6pt_mhv(twistor)
        H = H_res[0] if isinstance(H_res, tuple) else H_res
        
        # KLT
        permuted = [1, 2, 3]
        perms = list(permutations(permuted))
        
        def get_pt(order):
            denom = QQ(1)
            for k in range(6):
                denom *= twistor.get_angle(order[k], order[(k+1)%6])
            if denom == 0: return None
            num = twistor.get_angle(0, 1)**4
            return num / denom

        M_KLT = QQ(0)
        for alpha in perms:
            order_L = [0] + list(alpha) + [4, 5]
            AL = get_pt(order_L)
            if AL is None: continue
            for beta in perms:
                order_R = [0] + list(beta) + [5, 4]
                AR = get_pt(order_R)
                if AR is None: continue
                S = klt_momentum_kernel_6pt(list(alpha), list(beta), twistor)
                if S is None: continue
                M_KLT += AL * S * AR
        
        if H != 0 and M_KLT != 0:
            ratio = M_KLT / H
            ratios_mc.append(ratio)
            log_message(f"MC Sample {i}: Ratio = {float(ratio):.6f}")

    if len(ratios_mc) > 0 and len(set(ratios_mc)) == 1:
        log_message("[SUCCESS] Moment Curve Ratio CONSTANT!")

