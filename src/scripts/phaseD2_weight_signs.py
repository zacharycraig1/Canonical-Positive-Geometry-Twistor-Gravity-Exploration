import os
import sys
from sage.all import *

# Add src to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket
from src.chy_oracle.matrix_tree import hodges_minor_matrix_tree

def generate_positive_kinematics(n=6):
    """
    Generates kinematics on the moment curve: Z(t) = (1, t, t^2, t^3).
    Ordered t1 < t2 < ... < t6.
    Then converts to spinors.
    """
    # Sample ordered t
    ts = sorted([RR(randint(1, 100)) + RR(random()) for _ in range(n)])
    
    # Build Momentum Twistors Z
    Zs = [vector(RR, [1, t, t**2, t**3]) for t in ts]
    
    # Convert to Spinors
    # We need a map Z -> (lambda, tilde_lambda)
    # This requires choosing a splitting of infinity or similar.
    # Standard map:
    # lambda_i = (1, t_i) ? No, that's for 2D.
    # For 4D moment curve, we usually get specific spinors.
    # But simpler: Use the 4D moment curve directly in terms of Z, 
    # then compute brackets <ij> and [ij] from Z.
    # <ij> = det(Zi, Zj, I) ? No.
    # We need full 4-component twistors. Moment curve is in P3.
    # <i j k l> > 0.
    
    # But Hodges uses spinors.
    # We need "positive geometry" in spinor space.
    # Usually this means we map Z -> lambda, tilde_lambda.
    # Or we use the fact that for Moment Curve, <ij> > 0 (for j-i odd?).
    
    # Let's use the `src/klt_oracle/...` or `src/chy_oracle/amplitude_spinor`?
    # amplitude_spinor takes (lambda, tilde_lambda).
    
    # Alternative: Generate spinors directly such that <ij> > 0?
    # It is known that if we choose lambda_i = (1, t_i) with t1<...<tn,
    # then <ij> = t_j - t_i > 0 for i<j.
    # And tilde_lambda?
    # For real kinematics, tilde_lambda = conjugate(lambda) (for Euclidean?).
    # For split signature (2,2), tilde_lambda are independent real.
    # To be "positive", we usually work in split signature.
    
    # Let's try:
    # lambda_i = (1, x_i)
    # tilde_lambda_i = (1, y_i)
    # with ordered x, y.
    
    lambdas = []
    tilde_lambdas = []
    
    xs = sorted([RR(randint(1, 100)) + RR(random()) for _ in range(n)])
    ys = sorted([RR(randint(1, 100)) + RR(random()) for _ in range(n)])
    # Actually, for <ij> > 0, we need simple ordering.
    # <ij> = 1*xj - xi*1 = xj - xi > 0 for j > i.
    
    for x in xs:
        lambdas.append(vector(RR, [1, x]))
    for y in ys:
        tilde_lambdas.append(vector(RR, [1, y]))
        
    return lambdas, tilde_lambdas

def check_signs():
    print("Checking signs of weights w_ij = [ij]/<ij> on positive kinematics...")
    
    n_samples = 100
    consistent = True
    
    # Matrix to store sign patterns: S[i][j]
    # We expect pattern to be uniform across samples.
    sign_pattern = {}
    
    for s in range(n_samples):
        lambdas, tilde_lambdas = generate_positive_kinematics()
        
        # Compute weights
        for i in range(6):
            for j in range(i+1, 6):
                ang = ang_bracket(lambdas[i], lambdas[j])
                sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
                
                # Check ang > 0
                if ang <= 0:
                    # Should be positive by construction
                    print(f"Warning: <{i}{j}> = {ang} <= 0")
                    
                w = sq / ang
                sign = 1 if w > 0 else -1
                
                if (i,j) not in sign_pattern:
                    sign_pattern[(i,j)] = sign
                else:
                    if sign_pattern[(i,j)] != sign:
                        print(f"Sign flip at {i},{j} in sample {s}!")
                        consistent = False
                        
    if consistent:
        print("Success! Weight signs are deterministic on positive region.")
        print("Sign Pattern (i,j): Sign")
        for k, v in sorted(sign_pattern.items()):
            print(f"{k}: {v}")
            
        # Also check Tree Sum sign
        print("\nChecking Tree Sum signs...")
        tree_signs = []
        for s in range(10):
            lambdas, tilde_lambdas = generate_positive_kinematics()
            val = hodges_minor_matrix_tree(lambdas, tilde_lambdas)
            tree_signs.append(1 if val > 0 else -1)
            
        print(f"Tree Sum Signs: {tree_signs}")
        if all(x == tree_signs[0] for x in tree_signs):
            print("Tree Sum has uniform sign!")
        else:
            print("Tree Sum sign varies.")
            
    else:
        print("Signs are not consistent.")

if __name__ == "__main__":
    check_signs()



