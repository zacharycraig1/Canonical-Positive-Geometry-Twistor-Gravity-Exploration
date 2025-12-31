
from sage.all import *
from itertools import permutations

load('src/spinor_sampling.sage')
load('src/kinematics_map.sage')
load('src/klt.sage')
load('src/hodges.sage')

def debug_kernel():
    print("Debugging KLT Kernel...")
    
    # 1. Sample
    try:
        lambdas, tilde_lambdas = sample_spinor_helicity_conserving(n=6, seed=42)
    except: return
    
    Z_list = []
    x_current = matrix(QQ, 2, 2, 0)
    for i in range(6):
        lam = lambdas[i]
        til = tilde_lambdas[i]
        mu_0 = x_current[0,0]*lam[0] + x_current[0,1]*lam[1]
        mu_1 = x_current[1,0]*lam[0] + x_current[1,1]*lam[1]
        Z_list.append(vector(QQ, [lam[0], lam[1], mu_0, mu_1]))
        p_matrix = matrix(QQ, 2, 2, [lam[0]*til[0], lam[0]*til[1], lam[1]*til[0], lam[1]*til[1]])
        x_current += p_matrix
        
    twistor = MomentumTwistor(n=6, Z=Z_list, check_domain=False)
    twistor._compute_brackets()
    
    # Define direct adapter
    def mandelstam_direct(tw, i, j):
        # Ignore twistor, use captured spinors
        ang = lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
        sq = tilde_lambdas[i][0]*tilde_lambdas[j][1] - tilde_lambdas[i][1]*tilde_lambdas[j][0]
        return ang * sq
        
    # 2. Check Kernel Values
    alpha = [1, 2, 3] # Permutation of {1,2,3}
    beta1 = [1, 2, 3]
    beta2 = [3, 2, 1]
    
    # Use direct adapter
    val1 = klt_momentum_kernel_6pt(alpha, beta1, twistor, mandelstam_direct)
    val2 = klt_momentum_kernel_6pt(alpha, beta2, twistor, mandelstam_direct)
    
    print(f"S[123|123] = {val1}")
    print(f"S[123|321] = {val2}")
    
    if val1 == val2:
        print("Kernel value identical for different beta! (Suspicious)")
    else:
        print("Kernel values differ.")
        
    # Check if they are zero
    if val1 == 0:
        print("Value is zero!")
        
    # Trace inside kernel?
    # We can't modify src/klt.sage easily to print without reloading.
    # Let's verify what `klt_momentum_kernel_6pt` does manually here.
    
    print("\nManual Trace for S[123|123]:")
    # alpha = [1, 2, 3]
    # i=0: term = s_{0,1} + sum(j<0)... = s_{0,1}
    # i=1: term = s_{0,2} + theta(1,2)*s_{1,2}
    # i=2: term = s_{0,3} + theta(1,3)*s_{1,3} + theta(2,3)*s_{2,3}
    
    s01 = mandelstam_direct(twistor, 0, 1)
    s02 = mandelstam_direct(twistor, 0, 2)
    s03 = mandelstam_direct(twistor, 0, 3)
    s12 = mandelstam_direct(twistor, 1, 2)
    s13 = mandelstam_direct(twistor, 1, 3)
    s23 = mandelstam_direct(twistor, 2, 3)
    
    print(f"s01={s01}, s02={s02}, s03={s03}")
    print(f"s12={s12}, s13={s13}, s23={s23}")
    
    # Beta = [1, 2, 3]
    # theta(1,2) = 0 (1 is before 2). Wait.
    # Definition in klt.sage:
    # theta_beta(a, b) = 1 if a appears AFTER b in beta.
    # beta = [1, 2, 3]. pos(1)=0, pos(2)=1.
    # theta(1, 2): pos(1)=0 < pos(2)=1 -> 0. Correct.
    # theta(2, 1): 1 > 0 -> 1.
    
    # i=1 (alpha[1]=2): j=0 (alpha[0]=1).
    # theta(1, 2) = 0.
    # term = s02 + 0.
    
    # i=2 (alpha[2]=3): j=0 (1), j=1 (2).
    # theta(1, 3) = 0.
    # theta(2, 3) = 0.
    # term = s03.
    
    # S = s01 * s02 * s03.
    
    print(f"Manual S[123|123] = {s01*s02*s03}")
    
    # Beta = [3, 2, 1]
    print("\nManual Trace for S[123|321]:")
    # beta = [3, 2, 1]
    # pos(1)=2, pos(2)=1, pos(3)=0
    
    # i=0 (1): s01.
    # i=1 (2): j=0 (1). theta(1, 2). pos(1)=2 > pos(2)=1 -> 1.
    # term = s02 + s12.
    # i=2 (3): j=0(1), j=1(2).
    # theta(1, 3): pos(1)=2 > pos(3)=0 -> 1.
    # theta(2, 3): pos(2)=1 > pos(3)=0 -> 1.
    # term = s03 + s13 + s23.
    
    manual2 = s01 * (s02 + s12) * (s03 + s13 + s23)
    print(f"Manual S[123|321] = {manual2}")
    
    if manual2 != s01*s02*s03:
        print("Manual trace confirms dependence on beta.")
    else:
        print("Manual trace shows no dependence! (Implies s12=0 etc?)")

if __name__ == "__main__":
    debug_kernel()

