from sage.all import *

def sigma_mhv(kinematics, theta, eta, chi):
    """
    Computes the MHV rational solution sigma_a for n particles.
    Based on Du-Teng-Wu eq (2.15)/(3.18).
    
    sigma_a = <a eta> <theta chi> / ( <a chi> <theta eta> )
    
    Args:
        kinematics: SpinorKinematics object containing lambdas
        theta, eta, chi: Arbitrary reference spinors (2-vectors)
        
    Returns:
        List of sigma values [sigma_1, ..., sigma_n]
    """
    sigmas = []
    
    # Precompute constant factor <theta chi> / <theta eta>
    # <u v> = u0*v1 - u1*v0
    def angle(u, v):
        return u[0]*v[1] - u[1]*v[0]
        
    theta_chi = angle(theta, chi)
    theta_eta = angle(theta, eta)
    
    if theta_eta == 0:
        raise ValueError("Reference spinors theta and eta are collinear (<theta eta> = 0)")
    if theta_chi == 0:
        raise ValueError("Reference spinors theta and chi are collinear (<theta chi> = 0)")
        
    factor = theta_chi / theta_eta
    
    for a in range(kinematics.n):
        lambda_a = kinematics.lambdas[a]
        
        a_eta = angle(lambda_a, eta)
        a_chi = angle(lambda_a, chi)
        
        if a_chi == 0:
            # If <a chi> is 0, sigma_a goes to infinity. 
            # This is valid in CP1 but problematic for arithmetic in this form.
            # We might need to handle this or ensure reference spinors are generic.
            raise ValueError(f"Particle {a} spinor is collinear with chi")
            
        sigma_val = (a_eta * factor) / a_chi
        sigmas.append(sigma_val)
        
    return sigmas

def check_scattering_equations(kinematics, sigmas):
    """
    Verifies the scattering equations E_a = sum_{b!=a} s_ab / (sigma_a - sigma_b) = 0.
    
    Args:
        kinematics: SpinorKinematics object
        sigmas: List of sigma values
        
    Returns:
        Max absolute error (should be 0 for exact arithmetic)
    """
    n = kinematics.n
    residuals = []
    
    for a in range(n):
        E_a = 0
        for b in range(n):
            if a == b:
                continue
            
            s_ab = kinematics.s(a, b)
            sig_diff = sigmas[a] - sigmas[b]
            
            if sig_diff == 0:
                raise ValueError(f"Collision in sigma values: sigma_{a} == sigma_{b}")
                
            E_a += s_ab / sig_diff
            
        residuals.append(E_a)
        
    return residuals


