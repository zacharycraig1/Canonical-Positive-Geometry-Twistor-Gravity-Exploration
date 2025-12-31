# Cursor Next Steps: Toward a *canonical-form positive-geometry pushforward* statement
*(starting from this repo’s current “submission_pack” state)*

## 0) What we have **right now** (baseline)
You already have:
- A verified algebraic bridge **CHY MHV ↔ Hodges reduced determinant ↔ weighted Laplacian minor** (see `src/chy_oracle/laplacian_bridge.py`, `src/scripts/phaseE3_*`).
- A verified **all-minors / rooted-forest** expansion for the relevant Laplacian minors (see `src/scripts/phaseF4_all_minors_forest_expansion.py` and `src/chy_oracle/forest_enumerator.py`).
- Newton-polytope tooling (see `src/scripts/phaseD3_newton_polytope.py`, `src/scripts/phaseF5_newton_polytopes.py`).

**Missing for “positive geometry = canonical form = amplitude”:**
- A *precise* positive geometry **(X, X_{\ge0})** in an explicit ambient space.
- An explicit **canonical form** Ω(X) with *only logarithmic singularities* on boundaries.
- A **pushforward map** φ such that `φ_*(Ω(X)) = Ω_phys` (the physical amplitude/scattering form), with residues matching factorization.

This document is a concrete, implementable plan to get there.

---

## 1) Choose your first “target theorem”
You need a *single, crisp statement* you can try to prove / falsify.

### Target Theorem (T1, minimal but real)
For fixed roots `R={a,b,c}` and MHV reference spinors `(x,y)`, there exists a positive geometry
\[
(X_{n,R}, (X_{n,R})_{\ge0})
\]
together with a rational top-form Ω\_{n,R} (canonical form) and a pushforward map
\[
\phi_{n,R}: X_{n,R} \to \mathbb{P}\big(\mathcal{K}_n\big)
\]
into a projectivized kinematic space (or a natural subspace) such that
\[
\phi_{n,R*}(\Omega_{n,R}) \;=\; \Omega^{\text{MHV grav}}_n
\]
and the residues of Ω\_{n,R} on each boundary match factorization.

**Implementation rule:** Start with `n=4,5,6` and make this statement **explicit** and testable.

---

## 2) Two viable routes (do both in parallel, but keep them separate)

### Route A (Most aligned with your current “forest polytope” assets): **Toric / Newton-polytope positive geometry**
Key idea: your determinant minor is literally a **(rooted) Kirchhoff / forest polynomial** whose **Newton polytope** you already compute. There is a known positive-geometry story for **toric varieties associated to polytopes**: the positive part of a toric variety is a positive geometry with a canonical *dlog* form on the dense torus.

**Goal of Route A:** produce a clean pushforward formula:
- Start from the **torus** (positive coordinates),
- embed into a projective space via monomials (polytope lattice points),
- obtain a rational form whose denominator is the **forest polynomial / determinant**,
- then compose with your evaluation map `w_{ij}=[ij]/<ij>` (and the `C_i` factors) to reach the physics object.

This route is mathematically crisp and testable for small n, even if it is not yet “kinematic space ABHY-style”.

### Route B (More “physics-native”): **Scattering-form / kinematic-space geometry**
Key idea: treat amplitudes as differential forms on (projectivized) kinematic space and try to identify a positive region (polytope-like or arrangement-like) whose canonical form matches gravity in the MHV sector.
This route is harder, but it’s ultimately what you want if the claim is “gravity amplitude is canonical form of a positive geometry in kinematic space.”

---

## 3) Route A: Concrete deliverables and scripts (start here)

### A.1 Create a new phase directory
Create:
- `src/posgeom/` (new python/Sage module folder)
- `src/scripts/phaseG1_build_forest_polytope_data.py`
- `src/scripts/phaseG2_toric_coords_and_monomial_map.py`
- `src/scripts/phaseG3_canonical_form_from_facets.py`
- `src/scripts/phaseG4_pushforward_numeric_checks.py`

Add a single top-level runner:
- `src/scripts/PHASE_G_REPORT.md` (auto-generated summary; mirror style of `PHASE_E_REPORT.md`)

### A.2 Standardize “the polynomial”
Define a single canonical polynomial object:
- `F_{n,R}(z)` = the forest polynomial that equals the all-minors Laplacian minor in *formal* edge variables `z_e` (or `z_{ij}`).

**Do this:**
- In `src/posgeom/forest_polytope.py` implement:
  - `forest_polynomial(n, R) -> Sage polynomial in z_{ij}`
  - `forest_exponents(n, R) -> list of 0/1 exponent vectors`
  - sanity checks: degree must be `n-3` and monomials must correspond exactly to 3-rooted spanning forests.

**Must match existing computation:** cross-check against `phaseF4_all_minors_forest_expansion.py` results.

### A.3 Newton polytope = forest polytope (already mostly done, but formalize)
In `phaseG1_build_forest_polytope_data.py`:
- Build `P_{n,R} = Newt(F_{n,R})`
- Save:
  - vertices, facets, dimension, f-vector, lattice basis of affine span
- Output to: `.dcp_cache/phaseG/P_{n,R}_*.json`

**Pass criteria:**
- For `n=4,5,6`, polytope builds quickly and is stable across runs.
- Facets satisfy “matroid-style” inequalities (at least in the unrooted/tree limit). Don’t over-claim; just record them.

### A.4 Turn the polytope into a toric positive geometry (core step)
Implement in `src/posgeom/toric.py`:

1) Compute an integer basis `B` for the affine lattice of `P_{n,R}`:
- Map exponent vectors `m ∈ Z^{E}` into `Z^d` coordinates `u = B(m - m0)`.

2) Define the **monomial map** from a d-dimensional torus:
- Torus coords: `t=(t_1,...,t_d)`
- Each exponent vector `u` gives monomial `t^u`
- Embedding map: `t ↦ [t^{u_1}: ... : t^{u_N}]` where `{u_i}` enumerate lattice points you choose.

**Important simplification for first experiments:**
- For `n=4,5`, only use **vertices** (not all lattice points).
- For `n=6`, start with vertices; only extend to lattice points if needed.

3) The canonical form on the torus is:
\[
\Omega_{\text{torus}} = d\log t_1 \wedge \cdots \wedge d\log t_d
\]

Your job is to compute its pushforward into the projective embedding coordinates (or at least to compute the induced rational form in a convenient gauge).

### A.5 Compute the canonical form from facets (independent cross-check)
In `src/posgeom/canonical_polytope.py`, implement the **vertex–facet formula** for simple polytopes:

For a simple polytope with facet equations `ℓ_a(Y) ≥ 0`,
\[
\Omega(Y) \sim \sum_{v \in \text{vertices}} \frac{\mathrm{sgn}(v)\; d^dY}{\prod_{a \in \mathcal{F}(v)} \ell_a(Y)}
\]
where `F(v)` are facets incident to `v`.

**Why this matters:** it gives you an *explicit canonical form* you can numerically evaluate and compare against the toric pushforward form (A.4). If they disagree, your toric map/normalization is wrong.

### A.6 The key “pushforward check” you can actually run
In `phaseG4_pushforward_numeric_checks.py` implement:

- Pick random positive torus points `t ∈ (R_{>0})^d`.
- Compute image `Y = φ(t)` in embedding coords.
- Evaluate:
  1) `Ω_torus(t)` (as a number after choosing a coordinate chart / pulling back a test vector)
  2) `Ω_polytope(Y)` via vertex–facet formula (same chart)
- Confirm equality up to a global constant sign.

**You do not need full symbolic forms**; numeric equality at many random points is enough to validate the geometric model.

**Pass criteria:** for `n=4` and `n=5`, you get a consistent constant ratio across ≥200 samples.

---

## 4) Route A → Physics: connect toric/forest geometry to the gravity MHV function

Once Route A’s canonical-form machinery is correct, connect to your physics objects.

### A.7 Define the “evaluation map” into spinor kinematics
You currently evaluate weights as:
- `w_{ij} = [ij] / <ij>`
- `C_i = <i x><i y>`

Your determinant-minor identity uses:
- `det(L^{(R)}) / ∏_{k∉R} C_k^2` (plus normalization and `<ab>^8`).

Define a *single function* in `src/posgeom/physics_map.py`:
- `eval_edge_vars_from_spinors(lambdas, tildes, x, y) -> dict z_{ij} ↦ rational`
- Decide whether `z_{ij} = w_{ij}`, or `z_{ij} = w_{ij}/(C_i C_j)`, etc.
  - **Pick one** and keep it consistent everywhere.

### A.8 What you must test for “canonical form plausibility” (physics side)
In `phaseG5_physics_residue_tests.py` (new):
- Pull the candidate canonical form (Route A) forward to the chosen spinor subspace.
- Test **logarithmic** behavior on physical boundaries:
  - `<ij> → 0` (collinear)
  - `s_S → 0` (factorization) using kinematics that preserves momentum conservation

**Hard rule:** do *not* deform `<ij>` alone in a way that breaks momentum conservation and then conclude anything about factorization poles.

Use:
- Momentum-twistor sampling or a momentum-conserving deformation family (build it once; reuse).

**Pass criteria:**
- Poles are at most simple on each tested boundary (or you can explain why gravity should have higher order and adjust the claim).
- Residues match the expected product of lower-point MHV amplitudes (your existing factorization checker can be extended).

---

## 5) Route B: Kinematic-space scattering-form plan (in parallel, but later)
This is the “ultimate” claim if you want a gravity *kinematic-space* positive geometry.

### B.1 Build the kinematic scattering form object (n=6 first)
Implement a `src/posgeom/scattering_form.py` that can:
- enumerate cubic graphs (or use existing graph enumerator),
- attach propagators `s_I`,
- attach BCJ numerators (start with a known BCJ representation or numerically fit them),
- build the wedge-form sum.

Then:
- verify Jacobi in wedge space,
- verify pullback to planar subspace recovers YM partial amplitude.

### B.2 Try to get gravity from “square of YM” (experimental)
Song He and others have discussed gravity emerging from squaring structures of YM forms (e.g., “BCJ form → gravity”). Treat this as exploratory:
- compute residues of YM scattering form,
- square them in the sense appropriate for your representation,
- compare to your Hodges/Laplacian MHV amplitude.

This is not yet a canonical-form statement, but it can guide how a gravity geometry might be built.

---

## 6) Guardrails: when to stop / when to rewrite claims
### Stop conditions (do not proceed until fixed)
- If Route A fails already at `n=4` or `n=5`, pause and fix the toric/canonical implementation.
- If the physics residue tests show non-log poles generically, **do not** claim “positive geometry canonical form” yet; you must weaken the statement or change the geometry.

### Rewrite rules for paper (immediate)
Until Route A pushforward is validated and Route A→Physics residues pass:
- Replace “proves gravity amplitudes are positive geometries” with:
  - “identifies the Newton polytope / toric geometry naturally associated to the Hodges/Laplacian forest polynomial”
  - “investigates whether a canonical-form pushforward reproduces the physical amplitude”

---

## 7) File checklist (what Cursor should produce)
**New code**
- `src/posgeom/forest_polytope.py`
- `src/posgeom/toric.py`
- `src/posgeom/canonical_polytope.py`
- `src/posgeom/physics_map.py`
- `src/scripts/phaseG1_build_forest_polytope_data.py`
- `src/scripts/phaseG2_toric_coords_and_monomial_map.py`
- `src/scripts/phaseG3_canonical_form_from_facets.py`
- `src/scripts/phaseG4_pushforward_numeric_checks.py`
- `src/scripts/PHASE_G_REPORT.md` (auto-generated)

**Artifacts / caches**
- `.dcp_cache/phaseG/` with JSON dumps of polytopes + random-check logs

**Tests**
- `tests/test_phaseG_n4_pushforward.py`
- `tests/test_phaseG_n5_pushforward.py`

---

## 8) Minimal success definition for the next milestone
You can claim “we have a genuine canonical-form pushforward statement *in the toric/forest-polytope sense*” only if:

1) Route A: toric pushforward matches polytope canonical form numerically for `n=4,5` (constant ratio across many samples).
2) The induced form has only logarithmic singularities on toric boundaries (facet hyperplanes).
3) You can explicitly write down the map and the form in the paper (not just code).

Only **after** that do you attempt the full physics residue/factorization matching.

---

## 9) Suggested reading list to keep in the repo (add to `paper/refs.bib`)
Add bib entries for:
- Arkani-Hamed, Bai, Lam: “Positive Geometries and Canonical Forms” (canonical form axioms; toric examples)
- Arkani-Hamed et al.: “Scattering Forms and the Positive Geometry of Kinematics, Color and the Worldsheet” (kinematic-space forms)
- Chaiken: all-minors matrix-tree theorem (rooted forests enumerated by Laplacian minors)
- Any “Kirchhoff/forest polynomial ↔ matroid/graphic polytope” references you trust

(Keep this list in `CITATIONS.md` and mark clearly what is *prior art* vs what is *new*.)
