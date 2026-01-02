# Cursor Phase R: Turn “Forest Polytope” into an actual Pushforward / Canonical-Form Statement

**Context:** Phase P/Q built the *polytope side* (forest polytope inequalities, n=6 facet verification, canonical-form evaluator + tests).  
**What’s still missing:** (i) a **boundary dictionary** mapping facets ↔ physical limits, and (ii) an explicit **pushforward map** \(\Phi\) that makes the “canonical form = amplitude” claim precise.

This doc is a *Cursor-ready* checklist with concrete formulas, scripts to create, and pass/fail criteria.

---

## 0) Where we are now (do not redo)

### Proven / verified (n=6 baseline)
- Forest polytope for \(n=6, R=\{0,1,2\}\) has **108 vertices** and **22 facets**, matching the inequality model and expected counting.
- A canonical-form evaluator exists and passes simplex + square sanity tests.

### Known technical caveat (must address before publishing pushforward claims)
The canonical-form evaluator currently **sums `abs(term)`** to “force positivity/orientation consistency.” That’s ok as a *diagnostic* when \(W\) lies strictly in the dual-positive region, but it is **not** acceptable for residue/factorization statements or analytic continuation across poles.  
**Action:** replace `abs(term)` with a consistent sign choice (see Phase R2 below).

---

## 1) Explicit objects & formulas (to standardize everywhere)

### 1.1 Forest polynomial (combinatorial object)
For a root set \(R\subset [n]\), define
\[
F_{n,R}(z)\;=\;\sum_{F\in\mathcal{F}_R}\;\prod_{(i,j)\in E(F)} z_{ij},
\]
where \(\mathcal{F}_R\) is the set of spanning forests on \(K_n\) whose connected components are rooted at \(R\) (one root per component).

### 1.2 Weighted Laplacian & the all-minors identity
Define a weighted Laplacian \(L(z)\) by
\[
L_{ij} = -z_{ij}\quad (i\neq j),\qquad L_{ii}=\sum_{k\neq i} z_{ik}.
\]
Then (all-minors matrix-tree theorem):
\[
\det\big(L^{(R)}\big)=F_{n,R}(z),
\]
where \(L^{(R)}\) is the principal minor obtained by deleting rows/cols indexed by roots \(R\).

### 1.3 The physics map already in the repo
Using reference spinors \(x,y\) and spinors \((\lambda_i,\tilde\lambda_i)\),
\[
C_i=\langle i x\rangle\langle i y\rangle,\qquad 
w_{ij}=\frac{[ij]}{\langle ij\rangle},\qquad
z_{ij}=w_{ij}\,C_i C_j.
\]
This is the map \(\phi:\text{kinematics}\to \{z_{ij}\}\).

### 1.4 Forest polytope \(P_{n,R}\) (Newton polytope / incidence hull)
Let \(E\) be the edge set of \(K_n\), with variables \(x_e\) for each edge \(e\).  
The forest polytope is the convex hull of incidence vectors of rooted forests:
\[
P_{n,R}=\text{conv}\{\mathbf{1}_F \in \{0,1\}^{|E|}\;:\;F\in\mathcal{F}_R\}.
\]

A standard inequality description (via “super-root” spanning-tree reduction):
- \(0\le x_e\le 1\)
- \(\sum_{e\in E} x_e = n-|R|\)
- For each subset \(S\subset[n]\):
  - If \(S\cap R=\emptyset\): \(x(E(S))\le |S|-1\)
  - If \(S\cap R\neq\emptyset\): \(x(E(S))\le |S|-|S\cap R|\)
where \(x(E(S)):=\sum_{e\subset S} x_e\).

These are the families you must map to physics limits.

---

## 2) Phase R1: Build the **Boundary Dictionary** (facet → combinatorics → physics)

### Goal
Produce a machine-readable dictionary:
```json
facet_id -> {
  "family": "lower_bound" | "upper_bound" | "subset_rank",
  "edge": [i,j] (if bound),
  "subset_S": [list of vertices] (if subset_rank),
  "rhs": integer,
  "roots_in_S": integer,
  "predicted_limit": "...",
  "tests": {...}
}
```

### R1.1 Implement facet parsing (deterministic, no guessing)
Create:
- `src/pushforward/build_boundary_dictionary_n6.sage` (or `.py` if you already have Sage import conventions)

Algorithm sketch:
1. Load the **same edge ordering** used by `get_forest_polytope_inequalities(n, roots)`.
2. Build the Sage polyhedron `P` from the inequalities.
3. Extract *facet-defining inequalities* from `P.inequalities()` and normalize back to your “vec·x ≤ rhs” form.
4. For each facet, identify its family:
   - **Lower bound:** vec has one entry `-1` at edge e with rhs=0 ⇒ \(x_e\ge 0\).
   - **Upper bound:** vec has one entry `+1` at edge e with rhs=1 ⇒ \(x_e\le 1\).
   - **Subset rank:** vec entries are in `{0,1}` and correspond exactly to edges internal to some subset \(S\). Reconstruct \(S\) by:
     - Start with empty set.
     - For each edge coefficient 1, add its endpoints to a candidate set \(S\).
     - Verify that the coefficient pattern matches *exactly* the induced edge set on \(S\).
5. Store to `RESULTS/boundary_dictionary_n6.json`.

### R1.2 Attach “predicted physical meaning” **as a hypothesis**
Do **not** claim this is true yet; store hypotheses for testing:
- \(x_{ij}=0\) facet ↔ “drop edge” ↔ **\(z_{ij}\to 0\)** (e.g., \([ij]\to 0\) at fixed \(\langle ij\rangle\) and fixed \(C_i\)).
- \(x_{ij}=1\) facet ↔ “forced edge” ↔ **\(z_{ij}\to\infty\)** (often \(\langle ij\rangle\to 0\) at fixed \([ij]\), i.e. collinear pole).
- subset-rank facet ↔ “saturation of internal edges” ↔ candidate **multi-particle factorization channel** \(s_S\to 0\).

This is the *dictionary you will test next*.

---

## 3) Phase R2: Fix canonical-form orientation (remove `abs(term)`)

### Goal
Make `CanonicalFormEvaluator.eval_polytope(P,W)` return the **signed** canonical form.

### R2.1 Add a “reference orientation” for each simplex
Approach (pragmatic and works in practice):
1. Choose a reference point \(X_\star\) in the interior of the polytope (e.g., average of vertices).
2. For each simplex \(T=(Z_{i_0},...,Z_{i_d})\) in a triangulation:
   - Compute the oriented volume sign `sgn(det([Z_{i_k}-X_*]))`.
   - Multiply the simplex contribution by that sign.
3. Sum signed contributions (no abs).

### R2.2 Add a test that **fails** under abs but **passes** signed
Add:
- `src/tests/test_canonical_sign_consistency.py`

Test idea:
- Evaluate the form at a few \(W\) points where some terms are known to cancel (e.g., symmetric choices).  
- Ensure `abs`-version ≠ `signed`-version on at least one case.
- Ensure `signed` version is triangulation-independent.

---

## 4) Phase R3: Combinatorial factorization on facets (this is where “physics” can show up)

### Goal
For each subset-rank facet, compute the **face polynomial** and check if it factorizes into smaller forest polynomials (or a controlled sum of products).

### R3.1 Implement face-polynomial extraction
Create:
- `src/pushforward/face_polynomial_n6.sage`

For each facet inequality `vec·x ≤ rhs`:
1. Enumerate all forests \(F\in\mathcal{F}_R\) (you already have enumeration utilities).
2. For each forest, compute `dot(vec, incidence(F))`.
3. Keep forests saturating equality `dot = rhs`.
4. Form the face polynomial:
\[
F_{\text{face}}(z)=\sum_{F:\, \text{saturates facet}} \prod_{e\in F} z_e.
\]

### R3.2 Attempt symbolic factorization (Sage `factor`)
Compute `factor(F_face)` and match factors against:
- \(F_{|S|,R\cap S}\) on the induced complete graph on \(S\)
- \(F_{\text{contracted}}\) on the complement with \(S\) contracted to a node (this may require a “contraction polynomial” helper)

Output:
- `RESULTS/facet_factorizations_n6.json`
- A human-readable table in `docs/facet_factorization_summary.md`

**Important:** For MHV\(_6\), many physical channels have *zero residue* (you already saw this for a 3-particle channel), so “factorization” may show up as **vanishing** rather than a nontrivial product. Record both outcomes.

---

## 5) Phase R4: Connect facet factorization to **actual kinematic limits**

### Goal
Show that when you take a BCFW / soft / collinear limit in kinematics, the induced \(z_{ij}\) scaling selects the corresponding facet (or face) predicted by the boundary dictionary.

### R4.1 Implement a “limit probe” that records z-scalings
Create:
- `src/pushforward/kinematics_limit_probe_n6.sage`

For each physical deformation (already implemented in your harnesses):
- Soft limit: leg \(s\) soft (h=+2 and optionally h=-2).
- Collinear: \(\langle ij\rangle \to 0\) or \([ij]\to 0\).
- BCFW factorization: \(s_S\to 0\) for representative \(S\).

At each epsilon:
1. Compute spinors (already in harness).
2. Compute edge vars \(z_{ij}\) using the physics map.
3. Fit exponents \(\alpha_{ij}\) in \(z_{ij}\sim \epsilon^{\alpha_{ij}}\) by using multiple eps values.
4. Compare the pattern \(\alpha_{ij}\) to the facet families:
   - large negative exponent ≈ “goes to infinity”
   - positive exponent ≈ “goes to zero”
5. Emit a “best matching facet” score and store it.

Output:
- `RESULTS/limit_to_facet_matches_n6.json`

---

## 6) Phase R5: Only after R1–R4, attempt an actual **pushforward theorem** statement

### What counts as “publishable”
A publishable statement must pin down:
1. The positive geometry \(X_{\ge 0}\) (toric variety / polytope / hypersurface—be explicit).
2. The canonical form \(\Omega_X\) (as a differential form, not just a numeric evaluator).
3. The map \(\Phi\) from that geometry to kinematic space (or to the \(z\)-space you then map from kinematics).
4. A proof (or at minimum a precise conjecture + extensive exact checks) that:
\[
\Phi_*(\Omega_X)=\Omega_{\text{phys}}
\]
where \(\Omega_{\text{phys}}\) is a known scattering form / amplitude-associated form.

### Minimal next theorem target (concrete and checkable)
Start with a “two-step” theorem:
- (T1) **Combinatorial:** The forest polynomial equals the appropriate Laplacian minor (done).
- (T2) **Geometric:** The facets of \(P_{n,R}\) correspond to the subset-rank inequalities above and induce the observed face-polynomial factorization/vanishing.
- (T3) **Physical match:** Under the explicit kinematic map \(z_{ij}=([ij]/\langle ij\rangle) C_i C_j\), the kinematic limits (soft/collinear/BCFW) match those facet families.

This can be written as a serious “positive geometry dictionary” paper even before a full pushforward integral is nailed down.

---

## 7) “Don’t rediscover” guardrails (what is known vs what might be new)
- Matrix-tree / all-minors → rooted forest determinants: known.
- Spanning tree polytope inequality families: known.
- Hodges ↔ graphs/determinants in gravity: known.
Potentially new *for publication*: the *specific*, explicit **facet-by-facet dictionary** connecting the forest polytope’s face structure to **gravity’s limit behavior** and a precise **positive-geometry framing** that survives to \(n=6\) (and beyond).

---

## 8) Acceptance criteria (what to print in the final Phase R report)

Phase R is “done” when you can produce:
1. `boundary_dictionary_n6.json` with stable facet classification.
2. `facet_factorizations_n6.json` + readable summary.
3. `limit_to_facet_matches_n6.json` showing consistent matches for soft/collinear probes.
4. Canonical-form evaluator is **signed** and triangulation-independent.
5. A short `docs/pushforward_statement.md` that states exactly:
   - what \(X_{\ge 0}\) is,
   - what \(\Omega_X\) is,
   - what \(\Phi\) is,
   - what equality you claim, and in what sense (function vs form, pullback vs pushforward, domain restrictions).

---

## 9) Suggested command sequence (so Cursor can run it end-to-end)
From repo root:
1. `sage -python src/posgeom/compare_facets.py` (sanity: still 108 vertices / 22 facets)
2. `sage src/pushforward/build_boundary_dictionary_n6.sage`
3. `sage src/pushforward/face_polynomial_n6.sage`
4. `sage src/pushforward/kinematics_limit_probe_n6.sage`
5. `pytest -q` (ensure sign-consistency tests pass)

---

## 10) If you want to push beyond “reframing” into frontier territory
Once Phase R is solid, the high-upside directions are:
- Extend everything to **n=7** (first place truly new combinatorial complexity appears).
- Move beyond MHV: try to identify an NMHV numerator analogue and see if a related polytope/hypersurface structure emerges.
- Try to construct a genuine **scattering form** in kinematic space (ABHY-style), and show it matches your geometry under a map.

(Do these only after Phase R is clean; otherwise you’ll accumulate ambiguity.)
