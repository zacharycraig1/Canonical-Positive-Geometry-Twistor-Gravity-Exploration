import sys
from sage.all import *
import os
import random
import traceback
import json

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

from kinematics.spinors import SpinorKinematics
from chy.mhv_solution import sigma_mhv
from chy.scattering_eqs import detprime_phi
from chy.integrands import mhv_gravity_amplitude

# Load Hodges 
root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def collect_points_exact(n, mt_A, mt_B, theta, eta, chi,
                         target_pts=90, t_candidates=600, print_traces=3):
    """
    Stage A: Collect exact valid points with failure logging.
    """
    points = []
    fail = {"tilde_none": 0, "chy_exc": 0}

    t_list = [QQ(k) for k in range(1, t_candidates + 1)]
    # random.shuffle(t_list) # Deterministic for now is better for debug

    print(f"Starting collection of {target_pts} points...")
    
    for t_val in t_list:
        if len(points) >= target_pts:
            break

        # Construct Z(t)
        Z_t = [mt_A.Z[i] + t_val * mt_B.Z[i] for i in range(n)]
        mt_t = MomentumTwistor(n, Z=Z_t, check_domain=False)

        lambdas = [mt_t.get_lambda(i) for i in range(n)]
        tilde_lambdas = [mt_t.get_tilde_lambda(i) for i in range(n)]
        
        if any(tl is None for tl in tilde_lambdas):
            fail["tilde_none"] += 1
            continue

        k_t = SpinorKinematics(n, lambdas, tilde_lambdas)

        try:
            sigmas = sigma_mhv(k_t, theta, eta, chi)
            det_p = detprime_phi(sigmas, k_t)
            m_chy = mhv_gravity_amplitude(sigmas, k_t, det_p)
            points.append((t_val, m_chy))
            
            if len(points) % 10 == 0:
                sys.stdout.write(f"\rCollected {len(points)}/{target_pts}")
                sys.stdout.flush()
                
        except Exception:
            fail["chy_exc"] += 1
            if fail["chy_exc"] <= print_traces:
                print(f"\nFailure at t={t_val}")
                traceback.print_exc()
            continue

    print("")
    if len(points) < target_pts:
        raise RuntimeError(f"Only collected {len(points)}/{target_pts} points. Fail stats={fail}")

    print(f"Collected {len(points)} points. Fail stats={fail}")
    return points

def reconstruct_rational_adaptive(points, max_deg=80, holdout=20):
    """
    Stage B: Robust rational reconstruction with adaptive degree search.
    Optimized: Uses finite field probe to find degrees first.
    """
    print("Starting rational reconstruction...")
    
    # 1. Finite Field Probe
    # Use a large prime
    p_prime = 65521 
    GF_P = GF(p_prime)
    
    pts_gf = []
    for (ti, yi) in points:
        try:
            pts_gf.append((GF_P(ti), GF_P(yi)))
        except ValueError:
            # Skip points that blow up modulo p (denominator is 0 mod p)
            continue
            
    print(f"Probe: Mapped {len(pts_gf)} points to GF({p_prime})")
    
    # Run search in GF
    # ... (Reuse try_degrees logic but over GF)
    
    best_dp, best_dq = -1, -1
    
    # Need holdout for GF too
    valid_gf = pts_gf[:holdout]
    train_gf = pts_gf[holdout:]
    
    def try_degrees_gf(dp, dq):
        num_unknown = (dp + 1) + dq
        if len(train_gf) < num_unknown: return False
        
        rows, rhs = [], []
        for (ti, yi) in train_gf[:num_unknown]:
            row = []
            ti_pow = GF_P(1)
            for _ in range(dp + 1):
                row.append(ti_pow); ti_pow *= ti
            ti_pow = ti
            for _ in range(dq):
                row.append(-yi * ti_pow); ti_pow *= ti
            rows.append(row)
            rhs.append(yi)
            
        A = matrix(GF_P, rows)
        b = vector(GF_P, rhs)
        
        try:
            sol = A.solve_right(b)
        except Exception:
            return False
            
        # Validate
        p = sol[:dp+1]
        q = sol[dp+1:]
        
        def eval_rat(val):
            # manual eval
            num = sum(p[i] * val**i for i in range(len(p)))
            den = 1 + sum(q[j] * val**(j+1) for j in range(len(q)))
            if den == 0: return None
            return num / den
            
        for (ti, yi) in valid_gf + train_gf[num_unknown:]:
            val = eval_rat(ti)
            if val is None or val != yi:
                return False
        return True

    found = False
    for total in range(0, max_deg + 1):
        if total % 10 == 0: print(f"Probe scanning total degree {total}...")
        for dq in range(0, total + 1):
            dp = total - dq
            if try_degrees_gf(dp, dq):
                best_dp, best_dq = dp, dq
                found = True
                break
        if found: break
        
    if not found:
        raise RuntimeError("Probe failed to find rational function in GF.")
        
    print(f"Probe found degrees: deg(P)={best_dp}, deg(Q)={best_dq}")
    
    # 2. Exact Reconstruction
    print("Running exact reconstruction in QQ...")
    
    # Just call try_degrees for best_dp, best_dq over QQ
    # We need to reimplement try_degrees for QQ here or pull it out.
    # Inline for simplicity.
    
    R = PolynomialRing(QQ, 't')
    t = R.gen()
    
    train_pts = points[holdout:] # Use same split? Or all points? 
    # Use all points for max stability, but need holdout to verify.
    # Actually, we trusted the Probe degrees. 
    # Let's use holdout validation in QQ too to be sure.
    valid_pts = points[:holdout]
    
    dp, dq = best_dp, best_dq
    num_unknown = (dp + 1) + dq
    
    rows, rhs = [], []
    for (ti, yi) in train_pts[:num_unknown + 5]: # Use a few extra rows for overdetermined system
        row = []
        ti_pow = QQ(1)
        for _ in range(dp + 1):
            row.append(ti_pow); ti_pow *= ti
        ti_pow = ti
        for _ in range(dq):
            row.append(-yi * ti_pow); ti_pow *= ti
        rows.append(row)
        rhs.append(yi)
        
    A = matrix(QQ, rows)
    b = vector(QQ, rhs)
    
    try:
        # solve_right for overdetermined
        sol = A.solve_right(b)
    except Exception:
         # Try with exact count if overdetermined fails (singular)
         A = matrix(QQ, rows[:num_unknown])
         b = vector(QQ, rhs[:num_unknown])
         sol = A.solve_right(b)
         
    p = sol[:dp + 1]
    q_rest = sol[dp + 1:]

    P = R(sum(p[i] * t**i for i in range(dp + 1)))
    Q = R(1 + sum(q_rest[j] * t**(j + 1) for j in range(dq)))

    # Reduce
    g = gcd(P, Q)
    if g not in R: g = R(g)
    P = P // g
    Q = Q // g
    if Q.leading_coefficient() < 0:
        P = -P; Q = -Q

    return P, Q

def bracket_poly_from_twistors(mt_A, mt_B, i, j, t):
    # Build <ij>(t) in QQ[t] using lambda(t) = lambda_A + t*lambda_B.
    # Uses MomentumTwistor.get_lambda to stay consistent with your conventions.
    lamAi = vector(QQ, mt_A.get_lambda(i))
    lamBi = vector(QQ, mt_B.get_lambda(i))
    lamAj = vector(QQ, mt_A.get_lambda(j))
    lamBj = vector(QQ, mt_B.get_lambda(j))

    li0 = lamAi[0] + t*lamBi[0]
    li1 = lamAi[1] + t*lamBi[1]
    lj0 = lamAj[0] + t*lamBj[0]
    lj1 = lamAj[1] + t*lamBj[1]
    return (li0*lj1 - li1*lj0)

def analyze_denominator_by_gcd(Q, mt_A, mt_B, n=6):
    """
    Stage C: Pole analysis using exact polynomial division (GCD).
    """
    print("Analyzing denominator poles via GCD...")
    
    R = Q.parent()
    t = R.gen()

    powers = {}
    Qrem = Q

    bracket_polys = {}
    for i in range(n):
        for j in range(i + 1, n):
            aij = bracket_poly_from_twistors(mt_A, mt_B, i, j, t)
            if aij == 0:
                continue
            # Normalize to avoid scalar issues
            aij = aij / aij.leading_coefficient()  # monic-ish
            bracket_polys[(i, j)] = aij

    for pair, aij in bracket_polys.items():
        e = 0
        while True:
            # Check divisibility
            if Qrem % aij == 0:
                Qrem //= aij
                e += 1
            else:
                break
        if e > 0:
            powers[pair] = e

    print("Extracted <ij>(t) powers in Q(t):")
    sorted_powers = sorted(powers.items())
    for pair, e in sorted_powers:
        print(f"<{pair[0]}{pair[1]}>^{e}")

    if Qrem.degree() > 0:
        print("\nWARNING: leftover factor in Q(t) not explained by any <ij>(t):")
        print(Qrem.factor())
    else:
        print("\nAll denominator factors accounted for by <ij>(t).")

    return powers, Qrem

def run_harness():
    n = 6
    seed_A = 101
    seed_B = 102
    
    mt_A = MomentumTwistor(n, seed=seed_A)
    mt_B = MomentumTwistor(n, seed=seed_B)

    theta = vector(QQ, [1, 0])
    eta   = vector(QQ, [0, 1])
    chi   = vector(QQ, [1, 1])

    # Stage A
    try:
        points = collect_points_exact(
            n, mt_A, mt_B, theta, eta, chi,
            target_pts=100, t_candidates=600, print_traces=3
        )
    except RuntimeError as e:
        print(f"Stage A failed: {e}")
        return

    # Stage B
    try:
        P, Q = reconstruct_rational_adaptive(points, max_deg=100, holdout=30)
        print(f"P(t) degree: {P.degree()}")
        print(f"Q(t) degree: {Q.degree()}")
    except RuntimeError as e:
        print(f"Stage B failed: {e}")
        return

    # Stage C
    analyze_denominator_by_gcd(Q, mt_A, mt_B, n=n)

if __name__ == "__main__":
    run_harness()

