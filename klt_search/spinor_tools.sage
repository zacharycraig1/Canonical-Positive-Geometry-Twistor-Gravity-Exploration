
from sage.all import *

def spinors_from_twistors(twistor):
    n = twistor.n
    lambdas = []
    t_lambdas = []
    for i in range(n):
        lambdas.append(twistor.get_lambda(i))
        t_lambdas.append(twistor.get_tilde_lambda(i))
    return lambdas, t_lambdas

def twistors_from_spinors(lambdas, t_lambdas):
    n = len(lambdas)
    # We need to solve for mu_i given lambda_i and t_lambda_i.
    # Formula: t_lambda_i = ( <i,i+1> mu_{i-1} + <i+1,i-1> mu_i + <i-1,i> mu_{i+1} ) / (<i-1,i><i,i+1>)
    # Let D_i = <i-1,i><i,i+1>.
    # D_i * t_lambda_i = <i,i+1> mu_{i-1} + <i+1,i-1> mu_i + <i-1,i> mu_{i+1}
    
    # This is a system of 2n linear equations for 2n variables (components of mu).
    # Variables: mu_0^0, mu_0^1, ..., mu_{n-1}^0, mu_{n-1}^1.
    
    # Calculate coeffs
    def ang(i, j):
        return lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
        
    # Build Matrix
    # 2n x 2n
    M = matrix(QQ, 2*n, 2*n)
    RHS = vector(QQ, 2*n)
    
    for i in range(n):
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        
        D_i = ang(im1, i) * ang(i, ip1)
        rhs_val = D_i * t_lambdas[i]
        
        # Row 2*i corresponds to eq for t_lambda_i^0
        # Row 2*i+1 corresponds to eq for t_lambda_i^1
        
        RHS[2*i] = rhs_val[0]
        RHS[2*i+1] = rhs_val[1]
        
        # Coeffs for mu_{i-1}
        c_im1 = ang(i, ip1)
        M[2*i, 2*im1] = c_im1
        M[2*i+1, 2*im1+1] = c_im1
        
        # Coeffs for mu_i
        c_i = ang(ip1, im1)
        M[2*i, 2*i] = c_i
        M[2*i+1, 2*i+1] = c_i
        
        # Coeffs for mu_{i+1}
        c_ip1 = ang(im1, i)
        M[2*i, 2*ip1] = c_ip1
        M[2*i+1, 2*ip1+1] = c_ip1
        
    # Solve
    try:
        sol = M.solve_right(RHS)
    except:
        return None
        
    # Construct Z
    Z = []
    for i in range(n):
        z = vector(QQ, [lambdas[i][0], lambdas[i][1], sol[2*i], sol[2*i+1]])
        Z.append(z)
        
    return Z








