#!/usr/bin/env sage
from sage.all import *
load('src/hodges.sage')

def debug_hodges_scaling():
    print("Debugging Hodges Scaling")
    twistor = MomentumTwistor(n=6, check_domain=False)
    
    # Base
    res_base = hodges_6pt_mhv(twistor) # (val, reason)
    if res_base[0] is None: return
    val_base = res_base[0]
    
    # Scaled
    Z_scaled = [2 * z for z in twistor.Z]
    twistor_scaled = MomentumTwistor(n=6, Z=Z_scaled, check_domain=False)
    
    res_scaled = hodges_6pt_mhv(twistor_scaled)
    val_scaled = res_scaled[0]
    
    print(f"Total Ratio: {float(val_scaled/val_base):.4e}")
    # print(f"Power: {float(log(abs(val_scaled/val_base))/log(2)):.2f}")
    
    # Use explicit reference spinors
    # lx, ly are 2-vectors
    lx_base = vector(QQ, [1, 2])
    ly_base = vector(QQ, [3, 4])
    
    lx_scaled = 2 * lx_base
    ly_scaled = 2 * ly_base
    
    def get_components(tw, lx, ly):
        rows_delete = [0, 1, 2]
        cols_delete = [3, 4, 5]
        
        # We need to construct Phi manually as in hodges.sage
        n = 6
        Phi = matrix(QQ, n, n)
        
        lambdas = []
        for i in range(n):
            lambdas.append(vector(QQ, [tw.Z[i][0], tw.Z[i][1]]))
            
        for i in range(n):
            for j in range(n):
                if i != j:
                    ij_sq = tw.get_square(i, j)
                    ij_ang = tw.get_angle(i, j)
                    if ij_ang == 0: return None
                    Phi[i, j] = ij_sq / ij_ang
        
        # Diagonals
        def ang_vec(v1, v2):
             return v1[0]*v2[1] - v1[1]*v2[0]
             
        for i in range(n):
            diag_sum = 0
            ang_ix = ang_vec(lambdas[i], lx)
            ang_iy = ang_vec(lambdas[i], ly)
            for j in range(n):
                if i == j: continue
                ang_jx = ang_vec(lambdas[j], lx)
                ang_jy = ang_vec(lambdas[j], ly)
                
                term = Phi[i,j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
                diag_sum -= term
            Phi[i, i] = diag_sum
            
        rows_keep = [3,4,5]
        cols_keep = [0,1,2]
        Phi_minor = Phi[rows_keep, cols_keep]
        det = Phi_minor.det()
        
        def compute_c(indices):
            i, j, k = indices
            return 1 / (tw.get_angle(i,j)*tw.get_angle(j,k)*tw.get_angle(k,i))
            
        cI = compute_c(rows_delete)
        cJ = compute_c(cols_delete)
        
        return {
            "phi_elem": Phi[3,0],
            "det": det,
            "cI": cI,
            "cJ": cJ,
            "total": det * cI * cJ * (-1) # approx
        }

    comp_base = get_components(twistor, lx_base, ly_base)
    comp_scaled = get_components(twistor_scaled, lx_scaled, ly_scaled)
    
    if comp_base is None or comp_scaled is None:
        print("Failed to compute components")
        return

    for key in comp_base:
        b = comp_base[key]
        s = comp_scaled[key]
        if b == 0:
            print(f"{key}: Base is zero")
        else:
            ratio = float(s/b)
        pow = log(abs(ratio))/log(2)
        print(f"{key}: Power {float(pow):.2f}")

if __name__ == "__main__":
    debug_hodges_scaling()

