# CHY Line-Reconstruction Harness (Keep the Bot on Rails + Faster)

This note packages the recommended refactor of your “reconstruct amplitude on a line” script into a **3-stage, deterministic harness** that:

- **Stops the bot from getting lost** (no silent failures, no fragile interpolation).
- **Reconstructs `P(t)/Q(t)` robustly** (adaptive degree search + holdout validation).
- **Extracts pole structure without QQbar roots** (exact `gcd` / repeated division by `<ij>(t)` polynomials).
- Adds several **speed optimizations** and optional parallelism.

---

## Why the current script gets lost

### The 5 recurring failure modes
1) **Not enough valid points**: exceptions/invalid kinematics silently reduce sample count below what reconstruction needs.  
2) **Fragile kernel-based reconstruction**: taking `ker.basis()[0]` is unstable when kernel dimension > 1.  
3) **QQbar root hunting is slow & fragile**: pole attribution via algebraic roots is unnecessary.  
4) **Exceptions swallowed**: you lose information on where the pipeline is failing (σ collisions vs det′Φ singular vs integrand blow-ups).  
5) **Wrong bracket evaluation**: assuming `<ij>(t)` equals det of the first two twistor coordinates can be inconsistent with your actual `get_lambda` conventions.

---

## North-star workflow (do these in order)

### Stage A — Collect a fixed number of *valid* exact points
- Sample `t` values until you have exactly `N` valid points.
- Track failure counts by cause.
- Print the first few stack traces (or log them) so debugging is deterministic.

### Stage B — Robust rational reconstruction (adaptive degrees)
- Split into train/holdout sets.
- Search `(degP, degQ)` from small to larger until a candidate validates.
- Fix normalization (e.g. `q0 = 1`) so you don’t get scale ambiguity.

### Stage C — Pole structure without QQbar
- Build each quadratic `<ij>(t) ∈ QQ[t]` from **`get_lambda`** (not raw twistor coordinates).
- Factor `Q(t)` by repeated exact division / gcd.
- Any leftover factor indicates **non-`<ij>` poles** (i.e. spurious structure or wrong denominator hypothesis).

---

## Drop-in code (replace reconstruction + pole analysis)

> Keep your CHY point evaluation (`sigma_mhv`, `detprime_phi`, `mhv_gravity_amplitude`) as-is.  
> Replace interpolation + QQbar-root logic with the blocks below.

```python
from sage.all import *
import random, traceback

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
    return (li0*lj1 - li1*lj0).expand()

def collect_points_exact(n, mt_A, mt_B, theta, eta, chi,
                         target_pts=90, t_candidates=600, print_traces=3):
    # Collect exactly target_pts valid (t, m_chy) points.
    # Deterministic enough to debug: logs failures by category.
    points = []
    fail = {"tilde_none": 0, "chy_exc": 0}

    t_list = [QQ(k) for k in range(1, t_candidates + 1)]
    random.shuffle(t_list)

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
        except Exception:
            fail["chy_exc"] += 1
            if fail["chy_exc"] <= print_traces:
                traceback.print_exc()
            continue

    if len(points) < target_pts:
        raise RuntimeError(f"Only collected {len(points)}/{target_pts} points. Fail stats={fail}")

    print(f"Collected {len(points)} points. Fail stats={fail}")
    return points

def reconstruct_rational_adaptive(points, max_deg=60, holdout=20):
    # Find P/Q with minimal degrees that matches all points.
    # Fix q0=1 to remove scaling ambiguity; solve A x = b over QQ.
    R = PolynomialRing(QQ, 't')
    t = R.gen()

    pts = list(points)
    random.shuffle(pts)
    valid_pts = pts[:holdout]
    train_pts = pts[holdout:]

    def try_degrees(dp, dq):
        # Unknowns: p0..p_dp (dp+1), q1..q_dq (dq); q0 fixed to 1
        num_unknown = (dp + 1) + dq
        if len(train_pts) < num_unknown:
            return None

        rows, rhs = [], []
        # square system first (fast)
        for (ti, yi) in train_pts[:num_unknown]:
            row = []
            ti_pow = QQ(1)
            for _ in range(dp + 1):
                row.append(ti_pow); ti_pow *= ti
            ti_pow = ti
            for _ in range(dq):
                row.append(-yi * ti_pow); ti_pow *= ti

            rows.append(row)
            # Move + yi*q0 (with q0=1) to RHS:
            # P(ti) - yi*(1 + Σ qj ti^j)=0 -> P - yi Σ qj ti^j = yi
            rhs.append(yi)

        A = matrix(QQ, rows)
        b = vector(QQ, rhs)

        try:
            sol = A.solve_right(b)
        except Exception:
            return None

        p = sol[:dp + 1]
        q_rest = sol[dp + 1:]

        P = sum(p[i] * t**i for i in range(dp + 1))
        Q = 1 + sum(q_rest[j] * t**(j + 1) for j in range(dq))

        # Reduce
        g = gcd(P, Q)
        P //= g; Q //= g
        if Q.leading_coefficient() < 0:
            P = -P; Q = -Q

        f = P / Q

        # Validate on holdout + remaining training points
        for (ti, yi) in valid_pts + train_pts[num_unknown:]:
            if Q(ti) == 0:
                return None
            if f(ti) != yi:
                return None

        return (P, Q)

    # Search minimal total degree
    for total in range(0, max_deg + 1):
        for dq in range(0, total + 1):
            dp = total - dq
            out = try_degrees(dp, dq)
            if out is not None:
                P, Q = out
                print(f"Found rational function with deg(P)={P.degree()}, deg(Q)={Q.degree()}")
                return P, Q

    raise RuntimeError("No rational function found up to max_deg. Either not rational-in-t or degrees too high.")

def analyze_denominator_by_gcd(Q, mt_A, mt_B, n=6):
    # Exact pole attribution: factor Q(t) by quadratic <ij>(t) via repeated division.
    # No QQbar, no roots.
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
            aij = aij / aij.leading_coefficient()  # monic-ish
            bracket_polys[(i, j)] = aij

    for pair, aij in bracket_polys.items():
        e = 0
        while Qrem % aij == 0:
            Qrem //= aij
            e += 1
        if e > 0:
            powers[pair] = e

    print("Extracted <ij>(t) powers in Q(t):")
    for pair, e in sorted(powers.items()):
        print(f"<{pair[0]}{pair[1]}>^{e}")

    if Qrem.degree() > 0:
        print("\nWARNING: leftover factor in Q(t) not explained by any <ij>(t):")
        print(Qrem.factor())
    else:
        print("\nAll denominator factors accounted for by <ij>(t).")

    return powers, Qrem
```

### How to use it in your script

```python
def reconstruct_amplitude_on_line():
    n = 6
    mt_A = MomentumTwistor(n, seed=101)
    mt_B = MomentumTwistor(n, seed=102)

    theta = vector(QQ, [1, 0])
    eta   = vector(QQ, [0, 1])
    chi   = vector(QQ, [1, 1])

    points = collect_points_exact(
        n, mt_A, mt_B, theta, eta, chi,
        target_pts=90, t_candidates=600, print_traces=3
    )

    P, Q = reconstruct_rational_adaptive(points, max_deg=60, holdout=20)
    print("P(t) =", P)
    print("Q(t) =", Q)

    analyze_denominator_by_gcd(Q, mt_A, mt_B, n=n)
```

---

## Speed optimizations (big wins)

### 1) Avoid QQbar roots entirely (**major**)
The new approach uses **exact division** by `<ij>(t)` polynomials. This is faster and more reliable than `roots(QQbar)`.

### 2) Precompute bracket polynomials once (**major**)
`<ij>(t)` depends only on `(A,B)` and your lambda extraction; compute them once, reuse.

### 3) Tighten sampling: reject early (**medium**)
Before doing CHY:
- compute adjacent brackets quickly (`<i,i+1>(t)`), if any are zero → skip.
This prevents heavy CHY work at singular points.

### 4) Cache per-`t` evaluation outputs (**medium**)
Store:
- `sigmas`
- `det′Φ`
- `m_chy`
in a dictionary keyed by `t_val` so repeated runs reuse computed values (especially when debugging reconstruction degrees).

### 5) Parallel point collection (**optional big win, but be careful**)
Each `t` is independent. You can parallelize Stage A:
- `multiprocessing.Pool` can work, but Sage sometimes dislikes fork/serialization.
- prefer Sage’s `@parallel` if it’s stable in your environment.
- if parallelization is unstable, don’t use it—Stage A should be deterministic first.

### 6) Don’t build `MomentumTwistor(n, Z=...)` every time (**advanced**)
If `get_lambda` is linear in `Z` and `get_tilde_lambda` can be expressed with local minors, you can:
- precompute symbolic polynomials in `t` for each component of λ and \~λ,
- evaluate them quickly for each `t`,
- avoiding repeated object construction.
This can be a big speedup, but only after Stage A–C are stable.

---

## “Bot doesn’t get lost” checklist

1) **Never reconstruct** unless you have ≥ (unknowns + 20) valid points.  
2) **Never swallow exceptions** without counting them and printing the first ~3 traces.  
3) **Never use QQbar roots** for pole structure — use `gcd/division` against `<ij>(t)`.  
4) **Always validate** on holdout points not used to fit.  
5) **Cache failing seeds**: when any assertion fails, print `seedA, seedB, t_val` and write to a small JSON log.

---

## Interpretation guide

- If `Q(t)` factors entirely into powers of `<ij>(t)` polynomials:  
  your denominator hypothesis is consistent with “angle-bracket poles only” along that line.

- If leftover factors remain in `Q(t)` not divisible by any `<ij>(t)`:  
  either  
  1) your amplitude has additional poles along the line (spurious poles / wrong object), or  
  2) your `<ij>(t)` polynomials are not computed consistently with your kinematics extraction.

- If adaptive reconstruction never finds a rational function up to `max_deg`:  
  either the function is too high-degree or your pointwise CHY evaluation is inconsistent (most commonly σ choice or det′Φ normalization).

---

## Recommended next micro-goal for the agent

Before any “numerator geometry” claims, enforce:

**CHY pointwise oracle stability:**
- For 20 random kinematics points (not on a line), verify:
  - scattering equations residuals are exactly zero,
  - Möbius/reference changes don’t affect amplitude,
  - permutation tests pass.

Only then do line reconstruction / pole attribution.
