# Phase C3 “Northstar”: Twistor Boundary Geometry + Positivity (post–Phase C2)

You’ve already done the **right** thing in Phase C2:
- **All 3‑particle channels** at 6‑pt MHV show **finite scaling** (\u03b1≈0) under exact BCFW solving.
- That’s consistent with **helicity selection** (factorization residue vanishes), so **continuing to hunt 1/s poles here is not productive**.

## Should the agent continue?
**Yes — but pivot.**  
Phase C2 is now “done enough.” The agent should move from **channel-factorization diagnostics** to **twistor‑space singularities + the polynomial numerator geometry**.

The next objective is:

> **Characterize and constrain the polynomial numerator** (e.g. \u2115(Z), \u0124\u2115 = \u2115/\u27e801\u27e9^6) by its **boundary behavior** on the \u27e8ij\u27e9=0 divisors and by **symmetry/positivity**, and use that to build a geometry-native description (CHY ↔ canonical form language).

---

# Phase C3 Deliverables (what “done” looks like)

## D1 — Verified divisor/valuation table in twistor variables
A table of **pole orders of M** and **zero orders of \u2115, \u0124\u2115** on each boundary \u27e8ij\u27e9 → 0:
- cyclic neighbors: \u27e8i,i+1\u27e9
- non-neighbors: \u27e8i,j\u27e9
- special: \u27e801\u27e9

This must be done in a way that **does not confuse parameterization poles** with physical ones (i.e., control for spurious blowups).

## D2 — Minimal constraint set that uniquely fixes \u0124\u2115 up to scale (if possible)
Use linear constraints from:
- permutation symmetry (S6) or a known subgroup,
- boundary vanishing orders (from D1),
- soft/collinear limits *in spinor variables* (safe),
- CHY structural identities (e.g. reduced determinant relations).

## D3 — Positivity probe (geometry signal)
On positive kinematics (moment-curve twistors):
- sign/zero-locus behavior of \u0124\u2115,
- whether \u0124\u2115 defines a meaningful “boundary” (e.g. separates regions, has expected vanishing on facets),
- whether \u0124\u2115 factors or has structured irreducible pieces.

---

# Agent Instructions (action plan)

## 0) Freeze Phase C2 and refactor it into a library
**Goal:** stop re-running expensive experiments and prevent “agent gets lost in calculating.”

**Do:**
- Move the working pieces into `src/chy_oracle/`:
  - `amplitude_spinor.py` (your `hodges_6pt_mhv_spinor`)
  - `bcfw.py` (generic shift, solve `s_S(z)=0`)
  - `channels.py` (enumerate unique channels + recommended shifts)
  - `logging.py` (JSON lines + run metadata)
- Keep scripts in `src/scripts/` as thin wrappers.

**Success condition:** Every script is <200 LOC and imports from `src/chy_oracle/*`.

---

## 1) Build a robust valuation engine for twistor divisors \u27e8ij\u27e9 → 0
**Core idea:** instead of bisection, do **controlled epsilon deformations** and measure exponents.

### 1A) The deformation recipe (safe)
Given a seed kinematics (momentum twistors Z):
- pick a divisor \u27e8ij\u27e9,
- construct a 1-parameter family Z(\u03b5) that forces \u27e8ij\u27e9 ~ \u03b5 while keeping other brackets “generic.”

**Implementation suggestion (works well in practice):**
- Keep all Z_k fixed except Z_j \u2192 Z_j + \u03b5·V, with a random fixed vector V (QQ entries).
- Reject if this accidentally kills other brackets at leading order.

### 1B) Measure valuations
For each \u03b5 in a geometric sequence (e.g. 10^{-3}…10^{-10} in QQ or QQbar):
- compute:
  - M(\u03b5) using the CHY/Hodges oracle (spinor-only if possible),
  - D(\u03b5) (your discovered clearing denominator),
  - N(\u03b5) = M·D,
  - \u0124N(\u03b5)=N/\u27e801\u27e9^6 (when relevant).

Then estimate exponent k by log-slopes:
- k(M) = d log|M| / d log|\u03b5|
- and similarly for N, \u0124N.

**Deliverable:** a stable integer k for each divisor, with diagnostics.

### 1C) Script to create
Create: `src/scripts/phaseC3_valuations.py`

Outputs:
- `phaseC3_valuations.jsonl` with records:
  - seed, divisor, deformation details, observed exponents, rejection reasons

---

## 2) Re-verify the clearing denominator D(Z) using valuations (not interpolation)
You already have D(Z) from Phase B, but now you should “lock it” by checking:

- For each divisor:
  - k(M) + k(D) ≈ k(N) where k(N) is integer ≥ 0 (polynomial means no poles).
- Confirm the special cancellation order on \u27e801\u27e9:
  - N has factor \u27e801\u27e9^6 (per your Phase C report),
  - so \u0124N is finite/nonzero generically as \u27e801\u27e9→0.

This is faster and more reliable than full rational reconstruction on lines.

---

## 3) Symmetry: make it a *unit test*, not a narrative
Create: `src/scripts/phaseC3_symmetry_checks.py`

Tests (pass/fail with counters):
- permutation invariance of the **physical amplitude** M (up to sign convention)
- invariance of \u0124N under the symmetry subgroup you expect (start with the symmetries you’ve already seen, then expand)

**Important:** do symmetry tests in a parameterization that doesn’t break symmetry artificially.
- If momentum-twistor ordering introduces artifacts, test symmetry **in spinor-helicity** where possible.

---

## 4) Attempt reconstruction of \u0124N (degree 24) with constraints (optional but high value)
Full brute-force monomial bases explode. Instead:

### 4A) Use *constraint-driven* basis reduction
From valuations:
- if \u0124N must vanish on certain boundaries to given orders, enforce those as linear constraints.

### 4B) Use modular arithmetic to speed linear algebra
If your coefficients are rational and huge:
- sample kinematics over large finite fields (mod p),
- solve for coefficients mod p,
- CRT lift to QQ if needed.

This typically turns “impossible” exact solves into tractable ones.

### 4C) Candidate basis types
Try these *in order*:
1) products of 4-brackets only (degree 24 in Z-weight sense),
2) mixed 4-brackets × angle brackets, but grouped to match per-particle weights,
3) graph-generated ansatz (like you did before), but now with much smaller target degrees.

**Deliverable:** either (i) an explicit formula for \u0124N in a compact ansatz, or (ii) a proof that a given ansatz class cannot contain \u0124N.

---

## 5) Positivity probe (the actual “positive geometry” signal)
Create: `src/scripts/phaseC3_positivity_probe.py`

Run on:
- moment-curve positive twistors (your standard generator),
- random positive points inside the region.

Compute:
- sign of \u0124N,
- distribution of | \u0124N |,
- frequency of near-zeros under small perturbations,
- whether \u0124N seems to be a boundary-defining polynomial (changes sign across a codimension‑1 set).

**If it never changes sign** on the positive region, that’s itself meaningful.

---

# “Stop conditions” (to keep the agent from wandering)

The agent should **stop** Phase C3 once it has:
1) A stable valuation table for all \u27e8ij\u27e9 boundaries.
2) A symmetry test suite that passes on ≥ 200 random seeds.
3) A clear statement about \u0124N:
   - reconstructed formula **OR**
   - constrained family with dimension count + why it’s hard.
4) A short positivity report (plots optional, but statistics required).

After that, the next big step is “Phase D”: relate these divisor/positivity constraints to a **specific positive geometry model** (DCP charts / nested sets / canonical form matching).

---

# Commands (Docker / Sage)

From your repo root (PowerShell):

```powershell
docker run --rm -it -v ${PWD}:/workspace -w /workspace sage-cursor:latest `
  sage -python src/scripts/phaseC3_valuations.py
```

Interactive shell:

```powershell
docker run --rm -it -v ${PWD}:/workspace -w /workspace sage-cursor:latest bash
```
Then:
```bash
sage -python src/scripts/phaseC3_symmetry_checks.py
sage -python src/scripts/phaseC3_positivity_probe.py
```

---

# Notes for the bot (high-level guidance)

- **Phase C2 factorization** at 6‑pt MHV is now a solved diagnostic: it’s mostly helicity‑forbidden residues.
- The *real* geometry is in **twistor divisor structure** and the **polynomial numerator** (especially \u0124N).
- Avoid “blind fitting.” Prefer:
  - valuations (orders),
  - symmetry constraints,
  - soft/collinear behavior,
  - modular solves + lifting.

If the bot starts re-running long scans without adding new constraints, it’s off-track.
