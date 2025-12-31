from sage.all import *

class SpinorKinematics:
    """
    Manages spinor kinematics for n particles.
    Stores numerical values for lambda and tilde_lambda.
    Provides methods for angle <ij> and square [ij] brackets,
    and Mandelstam invariants s_ab.
    """
    def __init__(self, n, lambdas=None, tilde_lambdas=None):
        self.n = n
        self.lambdas = lambdas if lambdas is not None else []
        self.tilde_lambdas = tilde_lambdas if tilde_lambdas is not None else []
        
        if self.lambdas and self.tilde_lambdas:
             self._check_consistency()

    def _check_consistency(self):
        if len(self.lambdas) != self.n:
            raise ValueError(f"Expected {self.n} lambdas, got {len(self.lambdas)}")
        if len(self.tilde_lambdas) != self.n:
            raise ValueError(f"Expected {self.n} tilde_lambdas, got {len(self.tilde_lambdas)}")

    def angle(self, i, j):
        """Compute angle bracket <ij>."""
        # <ij> = det(lambda_i, lambda_j) = li_0 * lj_1 - li_1 * lj_0
        li = self.lambdas[i]
        lj = self.lambdas[j]
        return li[0] * lj[1] - li[1] * lj[0]

    def square(self, i, j):
        """Compute square bracket [ij]."""
        # [ij] = det(tilde_lambda_i, tilde_lambda_j)
        # Note: Convention varies, usually defined s.t. s_ij = <ij>[ji] = -<ij>[ij]
        # Or s_ij = <ij>[ij] if using different metric. 
        # Standard spinor helicity: s_ij = (k_i + k_j)^2 = 2 k_i.k_j = <ij>[ji].
        # Let's use [ij] = lti_0 * ltj_1 - lti_1 * ltj_0
        lti = self.tilde_lambdas[i]
        ltj = self.tilde_lambdas[j]
        return lti[0] * ltj[1] - lti[1] * ltj[0]

    def s(self, i, j):
        """Compute Mandelstam invariant s_ij = <ij>[ji]."""
        # Note the order [ji] for standard (+---) metric convention where s_ij = (k_i+k_j)^2
        return self.angle(i, j) * self.square(j, i)

    @classmethod
    def random_rational(cls, n, seed=None):
        """
        Generates random rational kinematics satisfying momentum conservation.
        Uses momentum twistors internally to guarantee conservation.
        """
        if seed is not None:
            set_random_seed(seed)
            
        # Generate random momentum twistors
        Z = []
        for _ in range(n):
            vec = vector(QQ, [QQ.random_element(num_bound=10, den_bound=10) for _ in range(4)])
            Z.append(vec)
            
        # Extract lambdas (first 2 components)
        lambdas = [vector(QQ, [z[0], z[1]]) for z in Z]
        
        # Extract tilde_lambdas using the formula derived from twistors
        # |i] = (<i,i+1> mu_{i-1} + <i+1,i-1> mu_i + <i-1,i> mu_{i+1}) / (<i-1,i><i,i+1>)
        # This guarantees sum |i><i| = 0
        
        def get_angle(i, j, Z_list):
            return Z_list[i][0]*Z_list[j][1] - Z_list[i][1]*Z_list[j][0]
            
        tilde_lambdas = []
        for i in range(n):
            im1 = (i - 1) % n
            ip1 = (i + 1) % n
            
            mu_i = vector(QQ, [Z[i][2], Z[i][3]])
            mu_im1 = vector(QQ, [Z[im1][2], Z[im1][3]])
            mu_ip1 = vector(QQ, [Z[ip1][2], Z[ip1][3]])
            
            ang_i_ip1 = get_angle(i, ip1, Z)
            ang_ip1_im1 = get_angle(ip1, im1, Z)
            ang_im1_i = get_angle(im1, i, Z)
            
            denom = ang_im1_i * ang_i_ip1
            if denom == 0:
                # Retry if we hit a singularity
                return cls.random_rational(n)
                
            num = mu_im1 * ang_i_ip1 + mu_i * ang_ip1_im1 + mu_ip1 * ang_im1_i
            tilde_lambdas.append(num / denom)
            
        return cls(n, lambdas, tilde_lambdas)

    def check_momentum_conservation(self):
        """Checks sum_i |i><i| = 0."""
        P_total = matrix(QQ, 2, 2)
        for i in range(self.n):
            # |i><i| = lambda_i * tilde_lambda_i^T
            # In index notation: (lambda_i)_alpha (tilde_lambda_i)_dot_alpha
            li = self.lambdas[i]
            lti = self.tilde_lambdas[i]
            term = matrix(QQ, [[li[0]*lti[0], li[0]*lti[1]], 
                               [li[1]*lti[0], li[1]*lti[1]]])
            P_total += term
        
        return P_total == 0


