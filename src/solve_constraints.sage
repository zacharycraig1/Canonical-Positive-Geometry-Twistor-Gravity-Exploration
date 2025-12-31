import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *
load("src/ansatz_builder.sage")
load("src/boundary_sampler.sage")
load("src/oracle_gravity_mhv6.sage")
load("src/dcp_eval.sage")

def build_and_solve():
    print("Phase 5: Solving for Unique Form")
    print("--------------------------------")
    
    # 1. Generate Orbits (Basis)
    orbits = generate_ansatz_orbits(verbose=True)
    orbit_keys = sorted(orbits.keys())
    
    print(f"\nBasis Dimension: {len(orbit_keys)}")
    
    # Map orbit index -> list of triples
    # Ansatz(c) = sum_k c_k * (sum_{trip in orbit_k} 1/(s_trip))
    
    # 2. Constraint Loop
    # We have 2 unknowns (c0, c1).
    # We need at least 2 independent constraints.
    # Constraints come from Residues on boundaries.
    # Boundary: s_S -> 0.
    # Residue of Ansatz = sum c_k * Res(Orbit_k)
    # Residue of Oracle = computed via boundary limit
    # Condition: sum c_k * Res(Orbit_k) = Res(Oracle)
    
    # We will simply evaluate Ansatz and Oracle at specific points NEAR the boundary
    # and fit: Ansatz(z) ~ Oracle(z) => Ansatz(z) * s(z) ~ Oracle(z) * s(z).
    
    # Constraint generation:
    # Pick a channel S.
    # Pick a point near s_S=0.
    # Equation: c0 * Val0 + c1 * Val1 = Target
    # Val_k = sum_{trip in orbit_k} 1/(s_trip) * s_S
    # Target = Oracle * s_S
    
    matrix_rows = []
    rhs_vec = []
    
    # Channels to probe
    # We should probe different types of channels to distinguish orbits.
    # Orbit 0 rep: ((1, 2), (1, 2, 3), (1, 2, 3, 4)) -> 2pt, 3pt, 2pt
    # Orbit 1 rep: ((1, 2), (1, 2, 3, 4), (1, 2, 5, 6)) -> 2pt, 2pt, 2pt (disjoint?)
    
    # Let's probe a 3-particle channel and a 2-particle channel.
    
    probe_channels = [
        (1, 2, 3), # 3-pt
        (1, 2),    # 2-pt
        (1, 3),    # another 2-pt
        (1, 2, 3, 4) # This is same as (5,6) 2-pt channel
    ]
    
    # Pre-generate channels list for evaluation
    C = generate_channels_n6()
    
    # Function to get orbit term value
    def eval_orbit(orbit_triples, l, lt):
        val = QQ(0)
        # We can construct a dummy candidate vector for dcp_eval
        # Or just sum manually since we have the triples
        # But we need to map spinor s_ij to channel s_S
        
        # Get channel values dict
        ch_vals = spinors_to_channels(l, lt)
        # Convert to 1-based canonical keys matching C format
        ch_vals_1 = {}
        full_set = set(range(1, 7))
        for k0, v in ch_vals.items():
            s1 = set([x+1 for x in k0])
            
            # Match generate_channels_n6 logic:
            # - Size 2: keep as is
            # - Size 3: canonicalize min(S, S^c)
            if len(s1) == 2:
                key = tuple(sorted(s1))
                ch_vals_1[key] = v
            elif len(s1) == 3:
                comp = full_set - s1
                key = min(tuple(sorted(s1)), tuple(sorted(comp)))
                ch_vals_1[key] = v
            else:
                # Should not happen for n=6 (spinors_to_channels returns size 2,3)
                pass
            
        for (i_idx, j_idx, k_idx) in orbit_triples:
            chA, chB, chC = i_idx, j_idx, k_idx # They are tuples
            
            # DEBUG: Print first lookup attempt
            if val == 0 and len(orbit_triples) > 0 and 'printed_debug' not in locals():
                print(f"DEBUG: Looking up {chA} in ch_vals_1")
                print(f"DEBUG: Keys in ch_vals_1 (sample): {list(ch_vals_1.keys())[:5]}")
                printed_debug = True
            
            sA = ch_vals_1.get(chA)
            sB = ch_vals_1.get(chB)
            sC = ch_vals_1.get(chC)
            
            if sA is None or sB is None or sC is None:
                # DEBUG missing keys
                if 'printed_missing' not in locals():
                    missing = []
                    if sA is None: missing.append(chA)
                    if sB is None: missing.append(chB)
                    if sC is None: missing.append(chC)
                    print(f"DEBUG: Missing keys in ch_vals_1: {missing}")
                    printed_missing = True
                continue
                
            if sA==0 or sB==0 or sC==0: continue
            
            val += 1 / (sA * sB * sC)
            
        return val

    # Sampling setup
    import random
    random.seed(int(555))
    
    epsilon = QQ(1)/1000
    
    print("\nGenerating constraints...")
    
    for i_ch, channel in enumerate(probe_channels):
        # Sample near boundary
        # Need base seed
        base_seed = random.randint(1, 10000)
        res = sample_spinor_helicity_conserving(n=6, seed=base_seed)
        l0, lt0 = res
        
        # Convert channel to set for boundary sampler
        if len(channel) > 3:
             # Canonicalize
             s = set(channel)
             comp = set(range(1, 7)) - s
             canon = min(tuple(sorted(s)), tuple(sorted(comp)))
             S_set = set(canon)
        else:
             S_set = set(channel)
             
        # Find shift
        # Need a pair (a, b)
        # a in S, b not in S
        # 0-based indices
        S_0 = set([x-1 for x in S_set])
        comp_0 = set(range(6)) - S_0
        if not comp_0: continue
        
        a = list(S_0)[0]
        b = list(comp_0)[0]
        
        z_pole = solve_boundary_shift(l0, lt0, S_0, (a, b))
        if z_pole is None: continue
        
        # Shift
        z = z_pole + epsilon
        l_eps, lt_eps = apply_shift(l0, lt0, a, b, z)
        
        # Compute s_S for this point
        # (We know it should be approx proportional to eps, but compute exactly)
        ch_vals = spinors_to_channels(l_eps, lt_eps)
        # Find value for S_set
        # Map S_set 1-based to 0-based key in ch_vals
        # Actually ch_vals has 0-based keys.
        # S_0 is 0-based.
        s_S_val = ch_vals.get(tuple(sorted(S_0)))
        if s_S_val is None:
             # Maybe key order diff
             continue
             
        # Evaluate Basis Functions (multiplied by s_S)
        row = []
        for o_idx, o_key in enumerate(orbit_keys):
            triples = orbits[o_key]
            val = eval_orbit(triples, l_eps, lt_eps)
            row.append(val * s_S_val)
            
        # Evaluate Oracle (multiplied by s_S)
        oracle_res = oracle_M6(l_eps, lt_eps)
        if oracle_res["value"] is None:
            print(f"Oracle failed at constraint {i_ch}")
            continue
            
        target = oracle_res["value"] * s_S_val
        
        print(f"Constraint {i_ch}: s_S={float(s_S_val):.2e}")
        print(f"  Row: {[float(x) for x in row]}")
        print(f"  Target: {float(target):.2e}")
        
        matrix_rows.append(row)
        rhs_vec.append(target)
        
    # Solve
    M = matrix(QQ, matrix_rows)
    b = vector(QQ, rhs_vec)
    
    print("\nLinear System:")
    print(M.str())
    print("RHS:", b)
    
    try:
        solution = M.solve_right(b)
        print("\nSOLUTION FOUND!")
        print(f"c0 (Orbit 0) = {solution[0]}")
        print(f"c1 (Orbit 1) = {solution[1]}")
        
        # Verify
        print("\nVerifying on generic point...")
        res = sample_spinor_helicity_conserving(n=6, seed=999)
        l_ver, lt_ver = res
        
        val_ansatz = QQ(0)
        for idx, sol_c in enumerate(solution):
             val_ansatz += sol_c * eval_orbit(orbits[orbit_keys[idx]], l_ver, lt_ver)
             
        val_oracle = oracle_M6(l_ver, lt_ver)["value"]
        
        print(f"Ansatz: {float(val_ansatz)}")
        print(f"Oracle: {float(val_oracle)}")
        print(f"Ratio: {val_ansatz/val_oracle}")
        
        if val_ansatz == val_oracle:
            print("EXACT MATCH!")
        else:
             print("Match up to scale?")
             
    except ValueError:
        print("\nNo solution or inconsistent system.")
    except Exception as e:
        print(f"\nSolver failed: {e}")

if __name__ == "__main__":
    build_and_solve()

