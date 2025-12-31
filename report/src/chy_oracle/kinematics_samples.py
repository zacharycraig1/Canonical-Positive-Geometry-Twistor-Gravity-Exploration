from sage.all import *
import numpy as np

class MomentumTwistor:
    """
    Momentum twistor representation for n particles.
    
    Stores Z_i ∈ QQ^4 for i = 0, ..., n-1.
    Precomputes all angle brackets <i j> and 4-brackets <i j k l>.
    """
    
    def __init__(self, n=6, seed=None, Z=None):
        self.n = n
        
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
                # Ensure non-zero vector
                while all(x == 0 for x in z):
                    z = vector(QQ, [
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11)),
                        QQ(np.random.randint(-10, 11))
                    ])
                self.Z.append(z)
        
        self._compute_brackets()
    
    def _compute_brackets(self):
        """Precompute all brackets once."""
        n = self.n
        
        # Angle brackets: <i j> = Z_i[0]*Z_j[1] - Z_i[1]*Z_j[0]
        self.angle = {}
        for i in range(n):
            for j in range(n):
                self.angle[(i, j)] = self.Z[i][0] * self.Z[j][1] - self.Z[i][1] * self.Z[j][0]
    
    def get_angle(self, i, j):
        """Get angle bracket <i j>."""
        return self.angle.get((i, j), QQ(0))
    
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

def sample_twistor(seed=None, n=6):
    """
    Creates a random MomentumTwistor instance.
    """
    return MomentumTwistor(n=n, seed=seed)

def sample_spinors_from_twistor(seed=None, n=6):
    """
    Generates momentum conserving spinors from a random twistor configuration.
    Returns (lambdas, tilde_lambdas).
    """
    twistor = MomentumTwistor(n=n, seed=seed)
    
    lambdas = [twistor.get_lambda(i) for i in range(n)]
    tilde_lambdas = []
    
    for i in range(n):
        lt = twistor.get_tilde_lambda(i)
        if lt is None:
            # If we hit a singularity (rare with random integers), try next seed
            if seed is not None:
                return sample_spinors_from_twistor(seed=seed+1, n=n)
            else:
                return sample_spinors_from_twistor(seed=None, n=n)
        tilde_lambdas.append(lt)
        
    return lambdas, tilde_lambdas

def assert_generic(lambdas):
    """
    Checks that no angle bracket <i,i+1> is zero (or other basic genericity).
    Returns True if generic, raises ValueError if not.
    """
    n = len(lambdas)
    for i in range(n):
        j = (i + 1) % n
        val = lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
        if val == 0:
            raise ValueError(f"Degenerate kinematics: <{i}{j}> = 0")
    return True

