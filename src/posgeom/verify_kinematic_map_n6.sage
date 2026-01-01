import sys
import os
import json
import random
from sage.all import QQ, vector

# Adjust path to import local modules
sys.path.append(os.getcwd())

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)

def get_s_map(u):
    s01, s02, s03, s04, s12, s13, s14, s23, s24 = u
    s = {}
    s[(0,1)] = s01; s[(1,0)] = s01
    s[(0,2)] = s02; s[(2,0)] = s02
    s[(0,3)] = s03; s[(3,0)] = s03
    s[(0,4)] = s04; s[(4,0)] = s04
    s[(1,2)] = s12; s[(2,1)] = s12
    s[(1,3)] = s13; s[(3,1)] = s13
    s[(1,4)] = s14; s[(4,1)] = s14
    s[(2,3)] = s23; s[(3,2)] = s23
    s[(2,4)] = s24; s[(4,2)] = s24
    
    s05 = -(s01 + s02 + s03 + s04)
    s[(0,5)] = s05; s[(5,0)] = s05
    s15 = -(s01 + s12 + s13 + s14)
    s[(1,5)] = s15; s[(5,1)] = s15
    s25 = -(s02 + s12 + s23 + s24)
    s[(2,5)] = s25; s[(5,2)] = s25
    
    C3 = -(s03 + s13 + s23)
    C4 = -(s04 + s14 + s24)
    C5 = -(s05 + s15 + s25)
    
    s34 = (C3 + C4 - C5) / 2
    s35 = (C3 - C4 + C5) / 2
    s45 = (-C3 + C4 + C5) / 2
    
    s[(3,4)] = s34; s[(4,3)] = s34
    s[(3,5)] = s35; s[(5,3)] = s35
    s[(4,5)] = s45; s[(5,4)] = s45
    return s

def get_s_subset(s_map, subset):
    val = 0
    sub = sorted(list(subset))
    for idx, i in enumerate(sub):
        for j in sub[idx+1:]:
            val += s_map[(i,j)]
    return val

def verify_map():
    print("Loading fit map...")
    fit = load_json("RESULTS/kinematic_map_n6_fit.json")
    facets = load_json("RESULTS/facet_dictionary_n6.json")
    
    X0 = vector(QQ, fit["X0"])
    # N is stored as list of lists (rows of N)
    N_rows = fit["N_rows"]
    # Reconstruct N as matrix. We need columns to be columns.
    # But sage matrix(rows) constructs from rows.
    # N_rows came from N.rows(). So matrix(N_rows) reconstructs N.
    from sage.all import matrix
    N = matrix(QQ, N_rows).transpose() # Wait. N.rows() gives rows of N.
    # In my fit script: N = N_ker.matrix().transpose().
    # So N is 15x11.
    # N.rows() gives 15 lists of length 11.
    # matrix(N_rows) gives 15x11 matrix. Correct.
    # Wait, in fit script I did: "N": json_friendly([list(row) for row in N.rows()])
    # So `N_rows` in JSON is list of 15 lists.
    # `matrix(QQ, fit["N_rows"])` creates 15x11 matrix. Correct.
    N = matrix(QQ, fit["N_rows"]) 
    
    t0 = vector(QQ, fit["t0"])
    T = matrix(QQ, fit["T"]) # 11x9
    
    matched_info = {m["facet_id"]: m for m in fit["matched_facets"]}
    
    print("Verifying on 50 fresh random points...")
    passes = 0
    fails = 0
    
    for i in range(50):
        u_vec = [QQ(random.randint(-20, 20)) for _ in range(9)]
        u = vector(QQ, u_vec)
        
        # Calculate t(s)
        # t = t0 + T*u
        t = t0 + T * u
        
        # Calculate X
        # X = X0 + N*t
        X = X0 + N * t
        
        # Check matched facets
        s_map = get_s_map(u_vec)
        
        point_pass = True
        
        for facet in facets:
            fid = facet["facet_id"]
            if fid not in matched_info:
                continue
            
            info = matched_info[fid]
            kappa = info["kappa"]
            subset = info["subset"]
            
            f_coeffs = [QQ(x) for x in facet["ineq_b_A"]]
            b_F = f_coeffs[0]
            A_F = vector(QQ, f_coeffs[1:])
            
            # L_F = b + A.X
            L_F = b_F + A_F.dot_product(X)
            
            # RHS = k * s_S
            s_S = get_s_subset(s_map, subset)
            rhs = kappa * s_S
            
            if abs(L_F - rhs) > 1e-10:
                print(f"Sample {i}: Facet {fid} Mismatch! L_F={float(L_F)}, RHS={float(rhs)}, Diff={float(L_F-rhs)}")
                point_pass = False
        
        if point_pass:
            passes += 1
        else:
            fails += 1
            
    print(f"\nVerification Results: {passes}/50 passed.")
    if fails == 0:
        print("Gate A Passed: Exact slack proportionality verified.")
    else:
        print("Gate A Failed.")

if __name__ == "__main__":
    verify_map()



