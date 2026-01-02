# Cursor Next Steps — Phase S (After Phase R/“K Output”)

**Goal:** turn the current *“forest polytope ↔ MHV gravity”* evidence into a **real canonical-form / pushforward theorem**, *without accidentally re‑deriving known results*, and with **full physical singularity coverage** (soft + all collinear + multi‑particle factorization).

---

## 0) Quick evaluation of the current output

You **did make progress**:

- You built a **facet (boundary) dictionary** for the `n=6, |R|=3` forest polytope and classified **22 facets** (with **3 still “unknown”**).  
- You replaced the “`abs(term)` patch” with a **signed triangulation rule** using a reference point `W_ref`, and you checked triangulation consistency in the toy regime.  
- You ran **limit probes** and found:
  - **Soft(leg 5)** and **collinear(3||4)** match **codimension‑1 facets**.
  - **collinear(4||5)** does **not** match any **codimension‑1 facet** in the current representation, suggesting it sits on a **higher‑codimension face** or requires a **different root/patch choice**.

➡️ This is *exactly* the right sort of bottleneck to hit: it’s where “toy positive geometry” becomes “physics geometry”.

**But you are not publish‑ready yet** because the **pushforward theorem requires residues/factorization on *all* physical boundaries**, and (4||5) is currently “missing” at codim‑1.

---

## 1) “Don’t rediscover known stuff” — the minimal related-work checkpoint

Before you write “new physics” claims, internalize the following:

### 1.1 Tree / forest expansions of Hodges are known
Hodges’ determinant for MHV gravity has long been interpreted via the **matrix‑tree theorem** and related tree/forest expansions.  
Key references to read and cite:
- Feng et al., **“Graphs, determinants and gravity amplitudes”** (2012).  
- Krasnov, **“Weighted Laplacians, cocycles and recursion relations”** (2013).

**So the novelty is *not* “gravity = sum over trees/forests.”**  
Your novelty (if it exists) is the **positive-geometry statement**: a *specific* polytope/canonical form plus a *specific* map whose pushforward gives the gravity amplitude.

### 1.2 Canonical forms for polytopes are standard positive-geometry tech
Arkani‑Hamed–Bai–Lam (2017) define positive geometries and canonical forms; polytopes are the basic example.  
Your work must go *beyond* “we can compute a canonical form”.

### 1.3 Matroids/trees/“amplitudes” are an active literature thread
There is substantial work connecting **oriented matroids**, canonical forms, and amplitude-like objects (including spanning-tree combinatorics) — notably by **Thomas Lam** and collaborators (2024–2025 era).  
This does **not** automatically subsume your claim, but it’s close enough that you must check it carefully.

**Action item:** add a “RELATED_WORK.md” section that explicitly distinguishes:
- known *tree/forest expansions* of MHV gravity,
- known *polytope canonical-form machinery*,
- what you claim is new: **which polytope + which map + which residue/physical interpretation**.

---

## 2) Phase S1 — Fix the physics probes (do limits on-shell + momentum-conserving)

Your current probing logic explicitly allows **off-shell** deformations (“we don’t strictly need momentum conservation”).  
That is *not reliable* for identifying physical boundaries.

### S1.1 Build momentum-conserving deformations
Implement **on-shell, momentum-conserving** families for limits:

**Soft leg `s`:**
- Scale spinors:  
  \( \lambda_s \to \epsilon\,\lambda_s,\ \tilde\lambda_s \to \epsilon\,\tilde\lambda_s \)  
  so \(p_s \sim \epsilon^2\).  
- Keep others fixed; verify total momentum conservation (re-solve if needed).

**Collinear `i || j`:**
Use a standard collinear parametrization:
- \(\lambda_i = \lambda_P\), \(\lambda_j = \lambda_P\)
- \(\tilde\lambda_i = z\,\tilde\lambda_P + \epsilon\,\tilde\eta\),  
  \(\tilde\lambda_j = (1-z)\,\tilde\lambda_P - \epsilon\,\tilde\eta\)  
with fixed \(z\in(0,1)\) and small \(\epsilon\).

**Multi-particle factorization** \(s_S\to0\):
Use a BCFW-style deformation or explicit solve for kinematics where
\( (\sum_{a\in S} p_a)^2 = \epsilon \).  
(For n=6 this can be done by constructing momentum twistors with one 4-bracket scaled.)

### S1.2 Re-run “facet triggering” with *correct* physics data
For each limit family:
- sample \(\epsilon\in\{10^{-3},10^{-4},10^{-5}\}\)
- compute the induced edge variables \(z_{ij}(\epsilon)\)
- define a score for each facet inequality \(A\cdot x \le b\) such as:
  - distance-to-hyperplane: \(\delta(\epsilon)=b-A\cdot x(\epsilon)\)
  - look for \(\delta(\epsilon)\sim \epsilon^{\alpha}\) with \(\alpha>0\)
- record the **set of facets** whose \(\delta\) goes to 0 fastest.

This replaces the earlier “exponents from ratios” hack and will correctly reveal **codim>1** behavior.

---

## 3) Phase S2 — Resolve the “missing collinear 4||5” singularity

Your report already suggests the right hypothesis: **4||5 is not a single facet** in the current patch.

### S2.1 Identify the *minimal face* hit by the 4||5 limit
Algorithm:
1. Run the on-shell collinear family for (4||5).
2. For each \(\epsilon\), compute \(x(\epsilon)\) (your polytope coordinates).
3. Evaluate all facet slack variables \(\delta_a(\epsilon)=b_a-A_a\cdot x(\epsilon)\).
4. Find facets with \(\delta_a(\epsilon)\to 0\).
5. The intersection of those facets is the candidate **face**.

Output:
- a JSON blob: `{ "limit": "4||5", "active_facets": [...], "face_dim_estimate": ... }`

### S2.2 Compute the *iterated residue* along the active facets
If the limit hits codim‑k face \(F\):
- compute \(\operatorname{Res}_{a_1}\cdots\operatorname{Res}_{a_k} \Omega_P\)
- compare to the expected physics factorization:
  - collinear splitting × lower-point amplitude

For gravity MHV, the collinear behavior is strong enough that you can check scaling and z-dependence numerically first, then exactify later.

---

## 4) Phase S3 — Patch dependence and root choices (symmetry restoration)

The polytope depends on the **choice of roots** \(R\).  
But the physical amplitude is **fully permutation symmetric**.

### S3.1 Repeat the polytope build for multiple root sets
For n=6, loop over all \(\binom{6}{3}=20\) root choices.
For each root set:
- build facet dictionary
- run the 6 collinear probes and 6 soft probes
- record which singularities appear as codim‑1 facets

Hypothesis to test:
- each root choice gives a *chart/patch* where only some singularities are codim‑1,
- the **union/gluing** over root choices recovers full singularity coverage.

### S3.2 Try a “symmetrized geometry” if patch-gluing is necessary
If no single root choice has all collinear channels as codim‑1:
- define \(\Omega_{\text{sym}} := \sum_{R} \Omega_{P_{n,R}}\)
- or define a single geometry built from the **common refinement** of the normal fans
- check whether the symmetrized canonical form matches the symmetric gravity amplitude.

This is the most plausible road to a *real* geometric statement.

---

## 5) Phase S4 — Make the canonical form engine bulletproof

### S4.1 Remove floats from classification and geometry
Facet classification currently stores `float(...)` copies.  
Replace with exact `QQ` storage everywhere; only convert to float for display.

### S4.2 Prove triangulation invariance with randomized triangulations
Add a test:
- choose 10 random triangulations of the same polytope (different interior points or reverse vertex orderings)
- evaluate \(\Omega_P(W)\) at 20 random dual points \(W\)
- assert the results are identical.

### S4.3 Guarantee the reference dual point is valid
You used `W_ref = (1,1,...,1)` heuristically.  
Add an explicit check:
- verify \(W_ref\cdot V > 0\) for all homogeneous vertices \(V\)  
so that no simplex denominator vanishes and sign calibration is meaningful.

---

## 6) Phase S5 — The *actual* pushforward equality test (n=6)

This is the core “physics” test.

### S5.1 Define the exact equality you are claiming
You need a single equation, with all variables and measures explicit:

> **Claim template:**  
> Let \(P\) be the rooted-forest polytope (specified).  
> Let \(\Omega_P\) be its canonical form on the appropriate dual/projective space.  
> Let \(\Phi\) map physical kinematics \((\lambda,\tilde\lambda)\) to that space (specified).  
> Then  
> \[
> \Phi^*(\Omega_P) = M_n^{\mathrm{MHV}}\,(\text{known prefactors})\, d^{\dim}\!\text{(kinematics)}.
> \]
> (Or: the *rational function* parts agree up to a constant.)

### S5.2 Implement the **ratio test** at random exact kinematics
For 50 random momentum-twistor samples (exact rationals, positive region):
1. compute physics amplitude \(M_6\) (Hodges / KLT oracle)
2. compute \(\Omega_P\) evaluated at \(W=\Phi(\text{kinematics})\)
3. check the ratio is constant.

**If the ratio is not constant:** the geometry/map is not yet correct (or you’re missing Jacobians).

---

## 7) Deliverables for the next agent run

Create a new folder `phase_s/` with:

- `phase_s/probe_limits_on_shell.sage`
- `phase_s/identify_active_facets.py`
- `phase_s/iterated_residue_check.py`
- `phase_s/root_sweep_n6.py`
- `phase_s/triangulation_invariance_test.py`
- `phase_s/ratio_test_pushforward_n6.sage`
- `phase_s/RESULTS.md` (auto-written by scripts with tables)

**Stop condition / success criterion:**  
You can explain *every* soft + collinear singularity as a **boundary (facet or face)** of a geometry, with residues matching factorization, and the pushforward ratio test is constant.

---

## 8) If you succeed — what is plausibly “new physics” here?

Not “trees in gravity” (known).  
Plausible novelty would be something like:

- **A specific positive geometry / gluing of charts** whose canonical form **pushes forward** to the full permutation-symmetric gravity MHV amplitude,
- a new **geometric explanation** of why Hodges/tree formulas have the exact residue/factorization structure,
- a bridge to **matroid/oriented-matroid canonical forms** that makes gravity’s structure look inevitable.

That’s frontier-adjacent. But you must earn it with the full singularity/residue story first.
