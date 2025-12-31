# CURSOR PROJECT UPDATE (v2) — Positive Geometry Gravity @ n=6 (OS^3 / DCP)
**Prepared:** 2025-12-29 (America/Los_Angeles)
**Owner:** Zach

## 0) WHAT THIS UPDATE DOES
This update turns the Cursor agent into a *targeted searcher* for the specific
next-phase objective implied by the Dec 25 checkpoint:

  A) Validate (or kill) the chart-level “HIT” observed in the leading-log proxy
     for the laminar family:
         N_seed = [(1,2), (1,2,3), (5,6)]   with boundary S=(1,2,3)
     by replacing the proxy pullback with a *faithful DCP chart pullback*.

  B) If strict constant-coefficient OS^3 invariants still fail globally, pivot
     immediately to a minimally-relaxed ansatz:
        “OS^3 invariants × kinematic-dependent symmetric numerators”
     and re-solve factorization + (optionally) Hodges oracle constraints.

This is designed to be copy/paste usable as a Cursor “agent mode” spec.

## 1) CURRENT FACTS (DO NOT RE-DISCOVER)
**Model:**
- n=6 channel arrangement {s_S=0 : 2 ≤ |S| ≤ 3} (15 channels total).
- OS^3 projective quotient basis dimension observed: ~2008.
- S6-invariant constant-coefficient subspace dimension observed: 2
    basis: w0, w1  (or your naming b0,b1).

**Empirical results to treat as “known”:**
- Channel-level strict L×R-only residue/support tests (v13–v15) found *no* hits
  across scanned 3-particle divisors.
- A chart-level leading-log proxy (v17b) produced a “HIT” for S=(1,2,3) on:
    N_seed = [(1,2),(1,2,3),(5,6)]
  with reported ratio a0/a1 = Infinity (i.e., coefficient of w1 forced to 0).
  This is only a *clue* until true DCP pullback is implemented.

**Interpretation to adopt:**
- “Pure constant OS^3 invariants imply strict factorization everywhere” is likely
  false; next phase must test true DCP pullback and/or allow numerators.

## 2) UPDATED SUCCESS CRITERIA (WHAT COUNTS AS PROGRESS)
**Tier-1 (Validation):**
- Implement a DCP chart pullback map for maximal nested sets N (or at least for
  charts extending N_seed) that produces explicit pulled-back dlog forms.
- Re-test the “HIT” chart family under the true pullback:
    If the HIT disappears → we have a strong no-go for the proxy.
    If it persists → escalate to multi-chart / multi-boundary validation.

**Tier-2 (Global factorization attempt):**
- For a candidate Ω (in whatever ansatz), for each 3-particle divisor s_S=0:
    Res_{s_S=0} Ω  ==  Ω_4^grav(L∪{P}) ∧ Ω_4^grav(R∪{P})
  in the strict L×R support sense *or* in a clearly-defined weakened sense
  (see §7) that still matches physical factorization.

**Tier-3 (Oracle / uniqueness):**
- On random MHV kinematic samples, candidate must match Hodges det'Φ up to an
  overall constant (or a small discrete convention class if that’s expected).
  Use this as a tie-breaker / uniqueness constraint.

## 3) IMMEDIATE TASK LIST (DO THESE IN ORDER)

### 3.1 Reproduce and freeze the proxy-HIT deterministically (fast)
- Run v17b (or the latest proxy code) and *hardcode the seed* N_seed:
    N_seed = [(12),(123),(56)]  (use your internal representation)
- Print and cache:
    - the pulled-back residue wedge terms
    - the derived L/R classification per u-variable
    - the exact a0,a1 solution (show that b=0 is forced)
- Save as:
    .dcp_cache/HIT_proxy_S123_seed12_56.sobj (or .json)

Why: once true pullback exists, you want an A/B test on the same chart seed.

### 3.2 Build a TARGETED nested-set sampler that starts from N_seed (critical)
Goal: stop wasting CPU on random nested sets; extend the known “good” laminar
seed to many maximal nested sets.

Implementation sketch:
- Input: building set G = connected flats (bitmasks), and required flats:
    Freq = {F(123)} plus optional {F(12), F(56)} if they exist in G.
- Build nested sets by incremental extension:
    N ← {F(123), F(12), F(56)} filtered to those actually in G
    while N not maximal:
        propose candidate flat H from G \ N that is compatible with N
        accept by nested-set axiom test (no forbidden antichain join in G)
        also enforce “coverage_ok” (avoid degenerate charts)
- Use heuristics:
    - prefer H that increases rank coverage (compute r(⋃N) quickly)
    - prefer H that is near-complementary to existing flats (reduce overlap)
    - keep 10–50 diverse completions per seed (beam search)

Deliverable:
- A function sample_nested_sets_seeded(required_flats, n_samples, rng_seed)
  returning a list of maximal nested sets N.

### 3.3 Implement TRUE(ish) DCP pullback (minimum viable version first)
You do NOT need the full wonderful model generality on day 1. You need a map
good enough to decide whether the proxy HIT is real.

Use this “adapted-basis + chain-rescaling” construction:

(1) Choose an “adapted basis” B of size r=10 for the dual space spanned by
    hyperplane forms {ℓ_e} (channels), adapted to N:
    - for each flat F in N ordered by inclusion, pick a subset B_F ⊂ F\(⋃children)
      of size Δr = r(F) - Σ r(children) (children = maximal proper flats in N under F)
    - ensure union(B_F) is independent and totals 10.
  (Greedy is fine; restart if fails.)

(2) Use {ℓ_b : b in B} as a basis for V*; then coordinates on V (kinematic space)
    can be taken as t_b = ℓ_b(p) for a point p ∈ V.

(3) Define DCP chart variables:
    - one “blowup” variable u_F for each flat F in N (or each irreducible flat),
    - plus “slope” variables v_b for b in B.

(4) Define the chart map (core idea):
    For each basis element b ∈ B, let Chain(b) = {F ∈ N : b ∈ F}, ordered by inclusion.
    Set:
        t_b = (∏_{F ∈ Chain(b)} u_F) * v_b
    Then for every channel e (not just b), define:
        ℓ_e(p) = Σ_{b∈B} c_{e,b} t_b   (precompute coefficients from linear algebra)
             = Σ c_{e,b} (∏_{F∈Chain(b)} u_F) v_b
    Factor out the *minimal* u-monomial common to terms with nonzero c_{e,b}:
        ℓ_e(p) = u_min(e) * (unit_e(u,v))   where unit_e is not divisible by u_min(e)

(5) Now dlog(ℓ_e) = dℓ_e/ℓ_e can be computed symbolically and expanded in
    dlog u_F and dv_b terms. Your OS^3 element is a linear combination of wedges
    of dlog(ℓ_e). After substitution, you can take residues along u_{F(123)}=0.

This is not a full formal proof of the DCP chart, but it is a *faithful blowup-like*
pullback that captures the “sum of monomials” feature missing in the proxy.

Priority: implement this for charts that extend N_seed and compare proxy vs true.

### 3.4 Re-run factorization tests on:
- the specific chart family completing N_seed (many completions),
- then a broader set of randomly sampled maximal nested sets containing F(123).

Report metrics:
- number of charts tested
- hits found (if any): (N, solution ratio a0:a1, residue term support summary)
- stability across completions and across other S

## 4) IF TRUE PULLBACK STILL YIELDS “NO GLOBAL HIT”: PIVOT PLAN (DO NOT STALL)

### 4.1 Minimal numerator ansatz (next lever)
You currently have invariant basis w0,w1 (constant coefficients). Replace:
    Ω = a0 w0 + a1 w1
with:
    Ω = f0(kin) w0 + f1(kin) w1
where f0,f1 are low-degree symmetric polynomials in channels.

Practical choices (start small):
- degree 0: constants (already done)
- degree 1: f_i = α_i * (Σ s_{ij}) + β_i * (Σ s_{123}-type)  [whatever independent sums exist]
- degree 2: use a small S6-invariant basis of quadratic polynomials in channel coordinates

Implementation:
- Build a small S6-invariant polynomial basis {p_k}.
- Let f0 = Σ x_k p_k,  f1 = Σ y_k p_k.
- Impose factorization constraints as linear equations on x_k,y_k by sampling many charts.
- Solve (over QQ if possible, else over several finite fields + CRT).

### 4.2 Add a Hodges/MHV oracle constraint if factorization underdetermines
If factorization leaves dimension >1, evaluate candidate on random MHV points and
require:
    candidate / Hodges  = constant
(or in a discrete finite set of allowed conventions if you’ve diagnosed that).

## 5) CODE ORGANIZATION (MAKE CURSOR AGENT FAST + SAFE)
Repo skeleton recommendation:
```
src/
  os/
    build_os_basis.sage
    s6_invariants.sage
  dcp/
    connected_flats.sage
    nested_set_sampler_seeded.sage   <-- NEW
    chart_pullback_true.sage         <-- NEW
    residue_factorization_tests.sage
  run/
    run_proxy_hit_repro.sage         <-- NEW
    run_true_pullback_scan_S123.sage <-- NEW
    run_global_scan_all_S.sage
.dcp_cache/
  flats_*.sobj
  closures_*.sobj
  hits_*.json
```

Hard requirements:
- Every expensive step caches to .dcp_cache with a versioned key.
- Every run prints progress every N trials (no silent stalls).
- Multiprocessing must be spawn-safe (if __name__ == "__main__": main()).

## 6) COMPUTE STRATEGY (AWS + LOCAL)
- Prefer RAM-heavy CPU boxes. If swapping occurs, results are meaningless.
- Keep .dcp_cache on local disk, not network mounts.
- Use tmux, htop, iotop. If disk-bound, reduce cache writes and increase RAM.
- Tune NPROC conservatively (often 8 is better than 32).

## 7) NOTE ON “STRICT L×R SUPPORT” VS PHYSICS
Your current strict test is a strong condition. If it continues to fail even for
good candidates, define and implement a slightly weaker but still-physical test:

**Option A (projective equivalence):**
- Allow residue to match L×R wedge up to addition of terms that are exact or
  vanish in the boundary quotient.

**Option B (sum-over-charts factorization):**
- Allow factorization only after summing residues over a compatible cover of
  charts (global object), not chart-by-chart purity.

If you adopt a weaker criterion, write it down formally and implement it as an
explicit linear projection, so “pass/fail” remains objective.

## 8) WHAT TO OUTPUT AFTER A RUN (STANDARDIZED REPORT)
Always save:
- git commit hash
- script version + config (TRIALS, NPROC, KEEP, MAX_CHARTS, seed requirements)
- charts tested, acceptance rates: cov_ok, nested_ok
- all HITs as JSON:
    {S, N, ratio, residue_summary, timestamp}
- if no hits, store a “no-hit certificate”:
    {S, n_charts, random_seed, constraints_used}

## 9) ONE-SHOT PROMPT FOR CURSOR AGENT MODE
Copy-paste this into Cursor:

“Goal: validate the proxy HIT for S=(1,2,3) on N_seed=[(12),(123),(56)] by
implementing a true DCP chart pullback (adapted basis + chain-rescaling) and
rescanning charts completing that seed. Do NOT waste time on unseeded random
nested sets until the seeded scan is done. Cache connected flats and closures,
print progress, and output HITs in JSON with nested-set data and ratios. If no
true-pullback HIT survives after ≥5000 seeded charts, immediately pivot to a
numerator ansatz f0(kin) w0 + f1(kin) w1 with low-degree S6-invariant polynomials
and re-solve factorization constraints.”
