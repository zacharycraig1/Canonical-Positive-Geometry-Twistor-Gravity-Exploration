
from sage.all import *
from itertools import permutations
import numpy as np

load('src/kinematics_map.sage')
load('src/klt.sage')

class MomentumTwistorAdapter:
    def __init__(self, Z):
        self.n = len(Z)
        self.Z = Z
        self._compute_brackets()
        
    def _compute_brackets(self):
        self.angle = {}
        # Simple angle bracket approximation if Z are twistors?
        # Wait, for KLT we need Mandelstam s_{ij}.
        # In Twistor space, s_{ij} = <i-1 i i+1 j> ... complicated?
        # No, s_{i,i+1} = <i-1 i i+1 i+2> / (<i-1 i><i+1 i+2>).
        # We need the full spinor map or just use Hodgest formula for s_ij?
        # s_{ij} = <i j> [i j].
        # We need to extract lambda, tilde_lambda from Z.
        # But for Positive Region, we sample Z directly.
        # Can we define s_{ij} from Z?
        # Yes, standard formula involves dual conformal invariants or infinity twistor.
        # We need an Infinity Twistor I.
        # Standard choice: I = (0,0,0,1) + (0,0,1,0)?
        # For "Positive Geometry", usually we work with Momentum Twistors and a fixed Infinity Twistor.
        pass

    def get_angle(self, i, j):
        # We need spinors lambda.
        # lambda_i is top 2 components of Z_i.
        L_i = self.Z[i][:2]
        L_j = self.Z[j][:2]
        return L_i[0]*L_j[1] - L_i[1]*L_j[0]

    def get_square(self, i, j):
        # [i j] depends on mu part and Infinity Twistor.
        # Standard convention:
        # [i j] = <i-1 i j-1 j> / (<i-1 i> <j-1 j>)
        # This requires 4-brackets.
        
        im1 = (i - 1) % self.n
        jm1 = (j - 1) % self.n
        
        num = self.get_four_bracket(im1, i, jm1, j)
        den = self.get_angle(im1, i) * self.get_angle(jm1, j)
        if den == 0: return 0
        return num / den

    def get_four_bracket(self, i, j, k, l):
        M = matrix(QQ, [self.Z[i], self.Z[j], self.Z[k], self.Z[l]])
        return M.det()

def mandelstam_from_twistors(tw, i, j):
    # s_ij = <i j> [i j]
    return tw.get_angle(i, j) * tw.get_square(i, j)

def sample_positive_grassmannian(n=6):
    # Sample Z such that all ordered minors <ijkl> > 0.
    # Use Matrix parametrization for Pos Grassmannian.
    # C = (1 M) where M is positive?
    # Or just sample random positive matrix and check minors?
    # Better: Parametrization.
    # For G(4,6), use cyclic connected minors?
    # Easier: Use the "Positive Matrix" generation method.
    # Start with positive matrix.
    
    while True:
        # Generate random integer matrix 4x6
        # To get positive minors, elements should vary?
        # Actually, standard way is specific parametrization variables.
        # Let's try brute force with heuristic:
        # Vandermonde-like matrix is positive?
        # Z_i = (1, t_i, t_i^2, t_i^3) with 0 < t_1 < ... < t_n.
        # This guarantees <ijkl> > 0 (Vandermonde determinant).
        
        ts = sorted([QQ(np.random.randint(1, 100)) + QQ(np.random.random()) for _ in range(n)])
        Z = []
        for t in ts:
            # Add some noise or scaling to make it generic
            # Z = [1, t, t^2, t^3]
            # Ideally we want generic positive.
            row = vector(QQ, [1, t, t**2, t**3])
            Z.append(row)
            
        # Check positivity
        adapter = MomentumTwistorAdapter(Z)
        all_pos = True
        for i in range(n):
            j = (i+1)%n
            k = (i+2)%n
            l = (i+3)%n
            if adapter.get_four_bracket(i,j,k,l) <= 0:
                all_pos = False
                break
        
        if all_pos:
            return adapter

def check_positive_region_signature(num_samples=10):
    print(f"Checking KLT Signature in Positive Region (G+(4,6))...")
    
    signatures = {}
    
    for k in range(num_samples):
        adapter = sample_positive_grassmannian(6)
        
        # Compute KLT Matrix
        permuted_set = [1, 2, 3]
        basis_perms = sorted(list(permutations(permuted_set)))
        
        matrix_data = []
        for alpha in basis_perms:
            row = []
            for beta in basis_perms:
                val = klt_momentum_kernel_6pt(list(alpha), list(beta), adapter, mandelstam_from_twistors)
                row.append(val)
            matrix_data.append(row)
            
        S = matrix(QQ, matrix_data)
        S_sym = (S + S.transpose()) / 2
        
        try:
            evals = S_sym.change_ring(RDF).eigenvalues()
            n_pos = sum(1 for e in evals if e > 1e-6)
            n_neg = sum(1 for e in evals if e < -1e-6)
            n_zero = sum(1 for e in evals if abs(e) <= 1e-6)
            sig = (n_pos, n_neg, n_zero)
            signatures[sig] = signatures.get(sig, 0) + 1
        except: pass
        
        if (k+1) % 5 == 0:
            print(f"Processed {k+1} samples...")
            
    print("\nSignatures in Positive Region:")
    for sig, count in signatures.items():
        print(f"  {sig}: {count}")

if __name__ == "__main__":
    check_positive_region_signature(20)








