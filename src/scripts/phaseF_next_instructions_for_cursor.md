# Phase F — From “Master Identity” to a Geometry-Native Proof (and n>6 Generalization)

**Status recap (what you have now):**
- You have an exact, numerically verified **Master Identity** expressing the 6-pt MHV gravity amplitude in terms of a **rank-3 weighted Laplacian** minor with CHY-consistent reference factors.
- You identified the **missing operation**: the corank-3 projection + universal prefactor that restores the physical double poles.
- You verified the **CHY ↔ Laplacian** bridge and the **positive spanning-tree geometry** signals.

**Phase F goal:** turn this into (i) a clean, citation-ready proof skeleton, and (ii) an implementation that *generalizes to n points* and exposes the **combinatorial positive-geometry object** (rooted forests / generalized permutohedra) that actually computes Hodges’ reduced determinant.

---

## 0) North Star (what “success” looks like)

By end of Phase F you should have:

1) A **proof-grade derivation** of the Master Identity from:
- CHY MHV solution → sigma-difference identity,
- CHY Jacobian → weighted Laplacian,
- Hodges reduced determinant → (n−3)x(n−3) minor due to corank 3,
- explicit normalization + reference-spinor cancellation.

2) A **combinatorial expansion** (no determinants) via the **All-Minors Matrix-Tree Theorem**:
- det of the (n−3)-minor = sum over **3-rooted spanning forests** of K_n,
- manifest positivity on the positive region (weights w_ij>0).

3) An **n=7 MHV verification** (and ideally n=8 spot checks):
- exact rational sampling,
- `hodges_npt_mhv_spinor` vs your Laplacian/forest formula,
- correct pole/limit behavior where residues *should* be nonzero.

---

## 1) Make the proof clean: isolate the “gauge transformations”

### Task F1.1 — Write the theorem statement and proof skeleton
Create:

- `src/notes/PHASE_F_THEOREM.md`

Include:
- precise definitions of `w_ij`, `C_i`, `L_tilde`, `Phi_tilde`,
- the corank-3 statement,
- the minor relation and prefactor,
- a short “why this is positive geometry” section (rooted forests + positivity of weights).

**Key identity to emphasize**
- `Phi_tilde = - D^{-1} L_tilde D^{-1}` with `D = diag(C_i)`.
- Therefore principal minors relate by explicit `prod(C_i)` factors.
- All reference dependence cancels in the final amplitude expression.

### Task F1.2 — Prove reference independence computationally (robustly)
Create:

- `src/scripts/phaseF1_reference_independence.py`

Algorithm:
1. Sample a valid (momentum-conserving) kinematic point (QQ).
2. Compute `M_phys` via `hodges_6pt_mhv_spinor`.
3. For ~50 random choices of reference spinors `(x,y)`:
   - build `L_tilde(x,y)`,
   - compute `M_rec(x,y)` via the Master Identity,
   - verify `M_rec(x,y)/M_phys == 1`.
4. Also check **gauge changes**:
   - rescale `x→αx`, `y→βy`, confirm invariance.

**Deliverable:** JSON log with worst deviation and any failure cases.

---

## 2) Replace determinants by combinatorics: All-Minors Matrix-Tree Theorem

This is the *real* “positive geometry object” behind your identity.

### Core fact (what you should implement)
Let `L` be a Laplacian with edge weights `a_{ij}` (symmetric).
For a root set `R` of size `k`, the principal minor `det(L_{V\R,V\R})` equals a sum over spanning **forests** with exactly `k` connected components, each containing exactly one root in `R`.

In your case:
- `n` particles,
- corank is 3,
- choose `R = {0,1,2}` (or any triple),
- the relevant minor size is `n−3`,
- so you sum over **3-rooted spanning forests** with `n−3` edges.

### Task F2.1 — Implement the rooted-forest expansion and match determinants
Create:

- `src/chy_oracle/forest_sum.py`

Functions:
- `k_rooted_forests_complete_graph(n, roots)`  
  returns an iterator of forests on K_n with components rooted at `roots` (k=3 fixed for now).
- `forest_sum_minor(lambdas, tildes, roots=(0,1,2), x=..., y=...)`  
  computes the forest sum with weighted edges.

Edge weight for your weighted Laplacian:
- `a_{ij} = w_ij * C_i * C_j` where `w_ij = [ij]/<ij>` and `C_i = <i x><i y>`.

Then verify:
- `det(L_tilde^{(roots)}) == forest_sum_minor(...)` exactly over QQ.

**Engineering note:** enumeration explodes for general n, but:
- for n=6, roots=3, forests have 3 edges → manageable,
- for n=7, forests have 4 edges → still manageable for correctness checks,
- for larger n, switch to determinant evaluation (fast) and keep forest enumeration only as validation.

### Task F2.2 — Extract the polytope: “rooted forest polytope”
The exponent vectors of the forest monomials live in R^{|E|}. For K_n the forest polytope is a generalized permutohedron variant.

Create:
- `src/scripts/phaseF2b_forest_newton_polytope.py`

Steps:
1. Represent each forest by its edge-incidence 0/1 vector in the |E| variables.
2. Build the polyhedron `P = Polyhedron(vertices=points)`.
3. Report: dim, #vertices, #facets.
4. Compare:
   - tree polytope (k=1 roots) vs 3-rooted forest polytope.

**Deliverable:** a small markdown report + cached `.sobj` of points.

---

## 3) Generalize to n points: the real breakthrough test

### Task F3.1 — Implement n-point MHV Hodges oracle (spinor-based)
If you don’t already have it, implement:

- `src/chy_oracle/amplitude_spinor.py`
  - `hodges_npt_mhv_spinor(lambdas, tildes, neg=(0,1), delete=(0,1,2))`

Make it:
- fully generic `n`,
- corank-3 reduced determinant computed by deleting any 3 rows/cols,
- normalized consistently.

### Task F3.2 — Implement the n-point weighted Laplacian reconstruction
Create:

- `src/chy_oracle/laplacian_bridge.py`
  - `weighted_laplacian(n, lambdas, tildes, x, y)`
  - `reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=(0,1,2), neg=(0,1))`

Then verify for:
- n=7 (first priority),
- n=8 (spot checks).

Create:
- `src/scripts/phaseF3_n7_verification.py`
- `src/scripts/phaseF3_n8_spotcheck.py`

**Pass criteria:** exact equality over QQ for ~50 seeds.

---

## 4) Pole structure & the “physical projection”: make it explicit

You already saw: the plain tree sum has the wrong pole order; the projected minor + prefactor fixes it.

### Task F4.1 — Automated valuation maps for n=7
Create:
- `src/scripts/phaseF4_n7_valuations.py`

For several divisors `<ij>→0`:
- measure slopes of:
  1) physical amplitude `M_phys`,
  2) weighted minor `det(L^{(roots)})`,
  3) full reconstructed expression (should match `M_phys`),
  4) optional: unweighted tree sum (for contrast).

**Goal:** identify which singularities are *intrinsic* to the forest geometry and which come purely from the universal prefactor.

### Task F4.2 — Factorization channels where residues should be nonzero (n=7+)
For 6-pt MHV many 3-particle channels had vanishing residues by helicity.
At n=7 you can find channels where both sides support ≥2 negative helicities (depending on split and internal helicity sums).

Create:
- `src/scripts/phaseF4b_n7_channel_factorization.py`

Use your existing BCFW engine:
- classify channels by predicted helicity-allowed residue vs forbidden residue,
- numerically confirm scaling.

---

## 5) Packaging: turn this into a paper-ready “result module”

### Task F5.1 — Consolidate all “final identities” into one place
Create:
- `src/notes/RESULTS_MASTER.md`

Include:
- The Master Identity,
- CHY sigma identity,
- Laplacian ↔ CHY Jacobian mapping,
- corank-3 argument,
- forest expansion statement,
- positivity region conditions (moment curve / ordered spinors).

### Task F5.2 — Minimal reproducibility harness
Create:
- `src/scripts/run_all_phaseF_checks.py`

It should:
- run n=6 sanity tests,
- run n=7 verification,
- run reference-independence checks,
- print a single “PASS/FAIL” summary.

---

## 6) Practical engineering notes (so Cursor stays fast and reliable)

- **Cache spanning trees/forests** as lists of edge-incidence vectors (not graph objects).
- Prefer **determinants** for n≥8; keep forest enumeration as correctness validation only.
- Keep everything in **QQ** where possible; only use floats for log-slope *display*.
- When slopes are needed, use exact fits on rationals:
  - compute `k = (log|f2| - log|f1|)/(log eps2 - log eps1)` with `RR` at the end.
- Always write JSON logs:
  - seed, parameters, failures, ratio, timings.

---

## 7) “What to do if something breaks”

- If `ValueError: reference spinor orthogonal...` → resample x,y (don’t skip the whole seed).
- If `tilde_lambda` reconstruction is fragile → for CHY tests, use your *momentum-conserving* sampler; avoid touching `<i,i+1>=0` charts unless the test is explicitly about that limit.
- If n=7 equality fails:
  1) check normalization (deleted rows/cols and norm factors),
  2) check helicity prefactor (which legs are negative),
  3) check sign conventions in `w_ij = [ij]/<ij>` and in the Laplacian diagonal.

---

## 8) The fastest “next 3 commits” (recommended order)

1) **Reference independence harness** (F1.2) + theorem write-up scaffold (F1.1).  
2) **All-minors forest sum** (F2.1) + determinant match at n=6.  
3) **n=7 verification** (F3.2) — if this passes, you’ve basically proven the general mechanism is correct.

---

### Quick checklist for Phase F “Done”
- [ ] `phaseF1_reference_independence.py` passes (many random x,y)
- [ ] `forest_sum_minor == det(L_minor)` at n=6 (exact)
- [ ] n=7: `M_laplacian == M_hodges` exact for multiple seeds
- [ ] valuation plots confirm pole mechanism persists
- [ ] `RESULTS_MASTER.md` contains a coherent proof narrative

