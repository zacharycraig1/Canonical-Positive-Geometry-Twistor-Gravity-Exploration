
import sys
from sage.all import *
load("correct_klt_proof.sage")

def generate_basis_graphs(target_degrees):
    # Optimized generator for specific degrees [8, 8, 2, 1, 1, 2]
    # We know 0-1 edge count is dominant.
    # Max possible <0,1> edges is 8.
    # Let's iterate k from 8 down to 0 for count of <0,1>.
    
    basis = []
    
    def solve(node_idx, current_degrees):
        if node_idx == 5:
            return [[]] if current_degrees[5] == target_degrees[5] else []
            
        needed = target_degrees[node_idx] - current_degrees[node_idx]
        if needed < 0: return []
        
        neighbors = range(node_idx + 1, 6)
        
        # Partition needed into neighbors
        import itertools
        
        # Optimization: for node 0, neighbor 1 is special.
        # We handle (0,1) edge count explicitly in outer loop if needed, or just let recursion handle it.
        # Let's stick to recursion but be clean.
        
        def partitions(n, k, limits):
            if k == 0:
                if n == 0: yield []
                return
            limit = limits[0]
            for i in range(min(n, limit) + 1):
                for rest in partitions(n - i, k - 1, limits[1:]):
                    yield [i] + rest

        limits = [target_degrees[j] - current_degrees[j] for j in neighbors]
        
        res = []
        for p in partitions(needed, len(neighbors), limits):
            edges = []
            next_deg = list(current_degrees)
            for i, count in enumerate(p):
                neigh = neighbors[i]
                next_deg[neigh] += count
                if count > 0:
                    edges.append(((node_idx, neigh), count))
            
            for s in solve(node_idx + 1, next_deg):
                res.append(edges + s)
        return res

    return solve(0, [0]*6)

def solve_modulo_prime():
    print("Solving Numerator coefficients modulo prime...")
    
    # 1. Basis
    target = [8, 8, 2, 1, 1, 2]
    basis = generate_basis_graphs(target)
    print(f"Basis size: {len(basis)}")
    
    # 2. Collect Data over GF(p)
    # We choose a large prime
    p = 1000000007 # 10^9 + 7
    F = GF(p)
    
    # We need to perform the hodges calculation in the field?
    # No, Hodges involves division (inverse).
    # We can compute in QQ and map to GF(p) at the end, 
    # OR compute entirely in GF(p).
    # Computing in GF(p) is faster and avoids expression swell.
    # But MomentumTwistor setup needs to be in GF(p).
    
    print(f"Working in GF({p})...")
    
    num_samples = len(basis) + 5
    X_rows = []
    y_vals = []
    
    count = 0
    seed = 1000
    
    while count < num_samples:
        # Generate random Z in GF(p)
        # We can't use 'sample_positive_Z_moment_curve' directly as it uses QQ and random floats
        # We implement a simple GF(p) sampler
        
        set_random_seed(seed + count)
        Z_entries = [[F.random_element() for _ in range(4)] for _ in range(6)]
        Z = [vector(F, z) for z in Z_entries]
        
        # Check domain (non-zero brackets)
        # Manually check required brackets
        valid = True
        
        # Check adjacent angles
        angles = {}
        for i in range(6):
            j = (i+1)%6
            ang = Z[i][0]*Z[j][1] - Z[i][1]*Z[j][0]
            if ang == 0: 
                valid = False; break
            angles[(i,j)] = ang
            
        if not valid: 
            seed += 1
            continue
            
        # We also need to run Hodges logic.
        # Can we run 'hodges_6pt_mhv' with Z in GF(p)?
        # Yes, if the function is generic.
        # But 'hodges_6pt_mhv' assumes QQ for some things?
        # It creates 'matrix(QQ, ...)' inside.
        # We need a version that accepts the field or infers it.
        # Let's monkey-patch or copy the logic briefly.
        
        # ... Copying core logic adapted for Field ...
        
        # Compute D
        D = F(1)
        for i in range(6): D *= angles[(i, (i+1)%6)]
        D = D**2
        
        # Compute M (Hodges)
        # Need Phi matrix
        Phi = matrix(F, 6, 6)
        
        # Reference spinors (random in F)
        lx = vector(F, [F.random_element(), F.random_element()])
        ly = vector(F, [F.random_element(), F.random_element()])
        
        # Check ref spinor validity
        # ... skip detailed check, if determinant fails we catch it
        
        # Helper for square brackets
        def get_sq(i, j):
            # [i j] = <i-1 i j-1 j> / (<i-1 i><j-1 j>)
            # 4-bracket
            im1, jm1 = (i-1)%6, (j-1)%6
            mat = matrix(F, [Z[im1], Z[i], Z[jm1], Z[j]])
            num = mat.det()
            den = angles[((i-1)%6, i)] * angles[((j-1)%6, j)]
            if den == 0: return None
            return num * den.inverse() # Field inverse
            
        fail = False
        for r in range(6):
            for c in range(6):
                if r != c:
                    sq = get_sq(r, c)
                    ang = Z[r][0]*Z[c][1] - Z[r][1]*Z[c][0]
                    if sq is None or ang == 0: 
                        fail = True; break
                    Phi[r,c] = sq * ang.inverse()
            if fail: break
        if fail: 
            seed += 1
            continue
            
        # Diagonal
        lambdas = [vector(F, [z[0], z[1]]) for z in Z]
        for r in range(6):
            ix = lambdas[r][0]*lx[1] - lambdas[r][1]*lx[0]
            iy = lambdas[r][0]*ly[1] - lambdas[r][1]*ly[0]
            if ix == 0 or iy == 0: 
                fail = True; break
            
            diag = F(0)
            for c in range(6):
                if r == c: continue
                jx = lambdas[c][0]*lx[1] - lambdas[c][1]*lx[0]
                jy = lambdas[c][0]*ly[1] - lambdas[c][1]*ly[0]
                term = Phi[r,c] * (jx*jy) * (ix*iy).inverse()
                diag -= term
            Phi[r,r] = diag
            
        if fail: 
            seed += 1
            continue
            
        # Reduced Det
        del_rows = [0, 4, 5] # Standard deletion
        rows = [1, 2, 3] # Kept
        
        Phi_red = Phi[rows, rows] # Symmetric
        det_red = Phi_red.det()
        
        # Norm factor
        # <0 4> <4 5> <5 0>
        # Note: 4,5 is adjacent. 0,4 and 5,0 are not.
        def get_a(u,v): return Z[u][0]*Z[v][1] - Z[u][1]*Z[v][0]
        
        n1 = get_a(del_rows[0], del_rows[1])
        n2 = get_a(del_rows[1], del_rows[2])
        n3 = get_a(del_rows[2], del_rows[0])
        norm = (n1*n2*n3)**2
        
        if norm == 0:
            seed += 1; continue
            
        M_val = det_red * norm.inverse()
        
        # Factor <0 1>^8
        M_val *= (get_a(0,1)**8)
        
        # N = M * D
        N_val = M_val * D
        
        # Monomials
        row = []
        for edges in basis:
            val = F(1)
            for (u,v), k in edges:
                val *= (get_a(u,v)**k)
            row.append(val)
            
        X_rows.append(row)
        y_vals.append(N_val)
        count += 1
        
    # Solve
    M_sys = matrix(F, X_rows)
    Y_sys = vector(F, y_vals)
    
    print(f"Solving {len(basis)}x{len(basis)} system in GF({p})...")
    
    # Check Rank
    rk = M_sys.rank()
    print(f"Rank: {rk}")
    
    if rk < len(basis):
        print("System is underdetermined (Basis is dependent). Finding one solution.")
        
    try:
        sol = M_sys.solve_right(Y_sys)
        print("Solution found!")
        
        # Print non-zero terms
        terms = []
        for i, c in enumerate(sol):
            c_int = int(c)
            # Map large field elements to negative integers
            if c_int > p // 2:
                c_int -= p
            
            if c_int != 0:
                # Format
                s = "".join([f"<{u}{v}>^{k}" if k>1 else f"<{u}{v}>" for (u,v),k in basis[i]])
                terms.append((c_int, s))
                
        terms.sort(key=lambda x: abs(x[0]), reverse=True)
        for c, s in terms:
            print(f"{c:+} * {s}")
            
    except Exception as e:
        print(f"No solution: {e}")

if __name__ == "__main__":
    solve_modulo_prime()







