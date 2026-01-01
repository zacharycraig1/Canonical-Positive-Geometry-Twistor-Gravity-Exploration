import sys
import os
import json
import random
import math
from sage.all import RR, QQ, vector, matrix

sys.path.append(os.getcwd())

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)

def compute_s_ij(lambdas, tildes, n=6):
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            li = lambdas[i]
            lj = lambdas[j]
            ti = tildes[i]
            tj = tildes[j]
            
            # <i j> = det(li, lj)
            bracket_angle = li[0]*lj[1] - li[1]*lj[0]
            # [j i] = det(tj, ti) -- usually [ij] is defined such that s = <ij>[ji]
            bracket_square = tj[0]*ti[1] - tj[1]*ti[0] 
            
            val = bracket_angle * bracket_square
            s[(i,j)] = val
            s[(j,i)] = val
    return s

def run_residue_check_fit():
    print("Loading geometry...")
    fit = load_json("RESULTS/kinematic_map_n6_fit.json")
    facets = load_json("RESULTS/facet_dictionary_n6.json")
    facet_dict = {f["facet_id"]: f for f in facets}
    
    X0 = vector(RR, fit["X0"])
    N = matrix(RR, fit["N_rows"])
    t0 = vector(RR, fit["t0"])
    T = matrix(RR, fit["T"])
    
    matched_info = {m["facet_id"]: m for m in fit["matched_facets"]}
    active_facets = [fid for fid, m in matched_info.items() if abs(m["kappa"]) > 1e-9]
    zero_facets = [fid for fid, m in matched_info.items() if abs(m["kappa"]) <= 1e-9]
    
    print(f"Active Facets (kappa != 0): {active_facets}")
    print(f"Zero Facets (kappa == 0): {len(zero_facets)} facets")
    
    # Probes
    probes = []
    
    # 1. Collinear (0,1) -> s_01 -> 0
    def probe_01(eps):
        ts = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        ts[1] = ts[0] + eps
        return ts
    probes.append(("Collinear 0,1 (s_01->0)", probe_01))
    
    # 2. Multi-particle (0,1,2) -> s_012 -> 0
    # s_012 = s01 + s02 + s12
    # Make p0, p1, p2 collinear? Or sum p0+p1+p2 -> 0 (soft)?
    # Collinear limit: p1 || p0, p2 || p0.
    def probe_012(eps):
        ts = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        ts[1] = ts[0] + eps
        ts[2] = ts[0] + 2*eps
        return ts
    probes.append(("Triple 0,1,2 (s_012->0)", probe_012))
    
    epsilons = [1e-2, 1e-4, 1e-6, 1e-8]
    n = 6
    
    for name, func in probes:
        print(f"\n--- Probe: {name} ---")
        
        last_slacks = {}
        
        for eps in epsilons:
            ts = func(eps)
            lambdas = {i: vector(RR, [1.0, ts[i]]) for i in range(n)}
            
            # Solve for conserved tildes
            # Random base tildes
            random.seed(42)
            tildes = {}
            for i in range(n-2):
                tildes[i] = vector(RR, [random.uniform(-1,1), random.uniform(-1,1)])
            
            L_mat = matrix(RR, [[lambdas[n-2][0], lambdas[n-1][0]], 
                                [lambdas[n-2][1], lambdas[n-1][1]]])
            rhs_vec_t0 = vector(RR, [-sum(lambdas[i][0] * tildes[i][0] for i in range(n-2)),
                                     -sum(lambdas[i][1] * tildes[i][0] for i in range(n-2))])
            rhs_vec_t1 = vector(RR, [-sum(lambdas[i][0] * tildes[i][1] for i in range(n-2)),
                                     -sum(lambdas[i][1] * tildes[i][1] for i in range(n-2))])
            try:
                sol_t0 = L_mat.solve_right(rhs_vec_t0)
                sol_t1 = L_mat.solve_right(rhs_vec_t1)
                tildes[n-2] = vector(RR, [sol_t0[0], sol_t1[0]])
                tildes[n-1] = vector(RR, [sol_t0[1], sol_t1[1]])
            except ValueError:
                print(f"Skipping eps={eps} due to singular spinor config")
                continue

            s_dict = compute_s_ij(lambdas, tildes, n)
            
            u_vals = [
                s_dict[(0,1)], s_dict[(0,2)], s_dict[(0,3)], s_dict[(0,4)],
                s_dict[(1,2)], s_dict[(1,3)], s_dict[(1,4)],
                s_dict[(2,3)], s_dict[(2,4)]
            ]
            u = vector(RR, u_vals)
            
            # Calculate X and Slacks
            X = X0 + N * (t0 + T * u)
            
            slacks = {}
            for fid, m in matched_info.items():
                f_coeffs = [float(x) for x in facet_dict[fid]["ineq_b_A"]]
                b_F = f_coeffs[0]
                A_F = vector(RR, f_coeffs[1:])
                L_F = b_F + A_F.dot_product(X)
                slacks[fid] = L_F
            
            # Identify minimal slack among active and zero sets
            min_act = min([abs(slacks[f]) for f in active_facets]) if active_facets else 0
            min_zero = min([abs(slacks[f]) for f in zero_facets]) if zero_facets else 0
            
            print(f"eps={eps:.1e} | MinActive={min_act:.2e} | MinZero={min_zero:.2e}")
            last_slacks[eps] = (min_act, min_zero)
            
        # Scaling Check
        if len(epsilons) >= 2:
            e1, e2 = epsilons[-2], epsilons[-1]
            s1, z1 = last_slacks[e1]
            s2, z2 = last_slacks[e2]
            
            if s1 > 1e-15 and s2 > 1e-15:
                slope = (math.log(s2) - math.log(s1)) / (math.log(e2) - math.log(e1))
                print(f"Scaling Active: ~ eps^{slope:.2f}")
            else:
                print("Scaling Active: Undefined (Zero/Noise)")
                
            if z1 > 1e-15 and z2 > 1e-15:
                slope_z = (math.log(z2) - math.log(z1)) / (math.log(e2) - math.log(e1))
                print(f"Scaling Zero: ~ eps^{slope_z:.2f}")
            else:
                print("Scaling Zero: Exact Zero (consistent with kappa=0)")

if __name__ == "__main__":
    run_residue_check_fit()
