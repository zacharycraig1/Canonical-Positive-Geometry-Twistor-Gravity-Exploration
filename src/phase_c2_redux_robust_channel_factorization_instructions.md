# Phase C2-Redux: Robust Channel Factorization Instructions (Avoid “Bisection to Infinity”)

**Problem observed:** Your bisection search for \(s_{123}(t)=0\) converged to a **pole** of \(s_{123}(t)\), where \(|s_{123}|\to\infty\), because \(s_{123}(t)\) is a **rational function** on a random twistor line and can change sign by passing through \(\infty\). This is a **chart/parameterization singularity**, not the physical channel.

This document gives **exact next steps** for the agent to implement immediately.

---

## North Star

> Replace “sign-bisection on \(s\)” with a method that targets **zeros of the numerator** of \(s\), while rejecting points where the denominator vanishes.

Use either:
- **Fix A (recommended on random lines):** Rationally reconstruct \(s_{123}(t)=P(t)/Q(t)\) and solve \(P(t)=0\) with \(Q(t)\neq0\).
- **Fix B (best overall):** Use a **BCFW-style shift** so \(s_{123}(z)\) is **linear** in \(z\) and solve for the root exactly.

---

## Important Physics Sanity Check (avoid false “failures”)

For **6-pt MHV gravity** with only two negative helicities, many multi-particle residues are **helicity-forbidden** at tree level (a factorization side may require a subamplitude with <2 minus helicities, which vanishes). Therefore it is plausible that at a true \(s_{123}\to0\) point:
- \(M\) is finite or vanishes,
- and \(s_{123}M\to 0\).

So the primary goal is:
1) **hit a true \(s_{123}=0\)** without chart poles,
2) measure the scaling of \(M\) vs \(s_{123}\),
3) and only then interpret whether a residue should be nonzero in this helicity sector.

---

## Fix A: Reconstruct \(s_{123}(t)=P(t)/Q(t)\) on a random twistor line

### Step A1 — Sample \(s_{123}(t)\) at rational points
- Choose a line: \(Z(t)=Z_A+t Z_B\) (random seeds OK).
- Sample \(t\in\{1,2,\dots,N\}\) (N ~ 40–80).
- For each point:
  - extract spinors \((\lambda,\tilde\lambda)\),
  - compute \(s_{123}(t)\) from momenta,
  - **skip** the point if extraction fails or if any adjacent \(\langle i,i+1\rangle=0\) (chart singularity).

Store a list of exact rational pairs:
\[
\{(t_i, s_{123}(t_i))\}.
\]

### Step A2 — Rational reconstruction (kernel solve)
Implement a kernel method like your earlier reconstruction harness:

```python
from sage.all import *

def rational_reconstruct_from_points(points, D):
    # points: list of (t_i, y_i) in QQ
    R = PolynomialRing(QQ, 't'); t = R.gen()

    rows = []
    for ti, yi in points[:2*D+2]:
        row = []
        # P coeffs
        ti_pow = QQ(1)
        for _ in range(D+1):
            row.append(ti_pow); ti_pow *= ti
        # Q coeffs
        ti_pow = QQ(1)
        for _ in range(D+1):
            row.append(-yi*ti_pow); ti_pow *= ti
        rows.append(row)

    M = matrix(QQ, rows)
    ker = M.right_kernel()
    if ker.dimension() == 0:
        raise RuntimeError("No kernel: increase D or collect more points")

    v = ker.basis()[0]
    p = v[:D+1]; q = v[D+1:]
    P = sum(p[i]*t**i for i in range(D+1))
    Q = sum(q[i]*t**i for i in range(D+1))

    # Optional: cancel gcd
    g = gcd(P, Q)
    if g.degree() >= 0:
        P //= g; Q //= g
    return P, Q
```

Choose \(D\) modestly (e.g. 10–25). Increase only if reconstruction fails.

### Step A3 — Solve \(P(t)=0\) while rejecting \(Q(t)=0\)
Compute roots of \(P\) in \(QQbar\); keep those where \(Q(r)\neq 0\):

```python
def safe_channel_roots(P, Q):
    roots = P.roots(ring=QQbar, multiplicity=False)
    good = []
    for r in roots:
        if Q(r) != 0:
            good.append(r)
    return good
```

### Step A4 — Validate “true zero” (not a pole) locally
For each candidate root \(t_*\):
- check nearby values \(t_*\pm\delta\) (small rational \(\delta\)) that:
  - \(|s_{123}|\) becomes small without blowing up,
  - adjacent brackets remain nonzero (avoid collinear/chart singularities),
  - the deformation does not force other invariants to vanish.

### Step A5 — Channel scaling test
Near \(t_*\), evaluate \(M\) and check scaling with \(s\):
- Fit \(M \sim s^\alpha\) from a few nearby samples.
- Record whether \(\alpha=-1\) (simple pole), \(\alpha\ge0\) (finite/vanishing), etc.

**Stop condition:** If you cannot find a safe root on a line (no real roots, or only roots where \(Q=0\)), resample a new line (new seeds).

---

## Fix B (preferred): BCFW shift so \(s_{123}(z)\) is linear and root is exact

This avoids “bisection-to-pole” completely.

### Step B1 — Extract a valid rational kinematics point
Start from one random twistor point \(Z\) that yields rational \((\lambda,\tilde\lambda)\).

### Step B2 — Apply a BCFW shift on legs \((a,b)\)
Use:
\[
\lambda_a \to \lambda_a + z\lambda_b,\qquad
\tilde\lambda_b \to \tilde\lambda_b - z\tilde\lambda_a.
\]
This preserves momentum conservation.

Implementation sketch:
```python
def bcfw_shift(lambdas, tilde_lambdas, a, b, z):
    lam = list(lambdas)
    tlam = list(tilde_lambdas)
    lam[a]  = lam[a]  + z*lam[b]
    tlam[b] = tlam[b] - z*tlam[a]
    return lam, tlam
```

### Step B3 — Solve for \(z_*\) with two evaluations
Compute \(s_{123}(z)\) at \(z=0\) and \(z=1\):
```python
s0 = s123_at_z(0)
s1 = s123_at_z(1)
c  = s1 - s0
z_star = -s0 / c
```
Now \(s_{123}(z_*)=0\) by construction (assuming linearity and \(c\neq0\)).

### Step B4 — Evaluate \(M(z)\) near \(z_*\) and infer scaling
Sample \(z=z_*+\delta\) for small rational \(\delta\), compute:
- \(s_{123}(z)\),
- \(M(z)\),
and fit \(\alpha\) in \(M\sim s^\alpha\).

---

## Pole Detector (keep this guard in all scripts)

Even after implementing Fix A/B, retain this guard:
- If \(|s|\) grows by > \(10^{k}\) while “searching,” abort interval (you are approaching a pole of \(s(t)\), not a zero).

**Never** accept a point where \(|s|\) becomes huge as “approaching the channel.”

---

## Required Logging (make debugging easy)

For every attempted channel test, log a JSON record including:
- seeds for \(Z_A,Z_B\) or base \(Z\),
- deformation type (line / BCFW),
- found root parameter (\(t_*\) or \(z_*\)),
- min adjacent bracket magnitude along samples,
- scaling exponent \(\alpha\),
- whether \(M\) diverged, finite, or vanished.

---

## Implementation Targets (files)

Create/update:

1. `src/scripts/phaseC2_factorization_robust.py`
   - Replace sign-bisection with Fix A or Fix B.
   - Add pole detector + logging.

2. `src/scripts/phaseC2_channel_scaling.py`
   - Given a root parameter, evaluate \(M\) vs \(s\) at several offsets and estimate \(\alpha\).

---

## Completion Criteria for This Patch

This patch is “done” when the agent can:
- reliably find a **true** \(s_{123}=0\) point without \(|s|\to\infty\),
- report the scaling exponent \(\alpha\),
- and show that results are stable under changing seeds/deformations.

Only after that should we interpret the outcome physically (helicity selection vs true factorization pole).

