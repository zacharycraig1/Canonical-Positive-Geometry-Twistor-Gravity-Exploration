# Cursor Next Steps (v3): Canonical-Form / Positive-Geometry Pushforward (Forest Polytope → MHV Gravity)

This note assumes the current Phase G artifacts:
- `forest_polytope.py` (forest polynomial + exponent list)
- `physics_map.py` (spinor-helicity → edge variables `z_ij`)
- `toric.py` (affine lattice basis + torus exponent coordinates)
- `canonical_polytope.py` (currently a placeholder / conceptual scratchpad)
- `PHASE_G_REPORT.md` (Phase G summary)

---

## 0) **Stop the novelty bleed** (what Phase G *actually* established)

### 0.1 What was verified
Phase G verified the identity
\[
F_{n,R}(\phi(\text{kinematics})) = \det(\tilde L^{(R)})
\]
for rooted forests and a weighted Laplacian, where `phi` is the edge-variable map `z_ij = ([ij]/<ij>) <ix><iy><jx><jy>`.
This is written explicitly in `PHASE_G_REPORT.md`.  

### 0.2 Why this is (almost certainly) not new
1. The equality “forest generating polynomial = Laplacian minor determinant” is the **all-minors matrix-tree theorem** (Chaiken and many later sources).
2. The fact that **MHV gravity** has a **spanning-tree / determinant** representation is also known (NSVW tree formula; Hodges determinant; and work explicitly connecting them via matrix-tree / graph determinants).

**Action:** rewrite any “new discovery” claims to: *we implemented and validated known combinatorial identities in our exact arithmetic / kinematic embedding pipeline*, and treat them as infrastructure for the new (still-missing) pushforward statement.

---

## 1) Patch the main gap: `canonical_polytope.py` is not doing the job

Right now `canonical_polytope.py` is explicitly unsure about what space the canonical form lives on and which pushforward is being checked.
Before writing more code, decide which of the following **two** pushforward statements you want.

---

## 2) Choose the target pushforward statement (pick one, and make it precise)

### Option A (toric pushforward → polytope canonical form; math-true, physics-unclear)
This is the statement in Arkani-Hamed–Bai–Lam (ABL) that the canonical form of a polytope (or its projective cone) can be obtained as a pushforward from a toric positive geometry / moment map.

- Pros: “book theorem” territory; you can implement and verify.
- Cons: does **not** by itself produce the gravity amplitude unless you add a *second* map tying the polytope form to kinematics.

### Option B (polytope in **kinematic** space; physics-true target)
This mirrors ABHY/associahedron logic:
- Build a polytope **directly in kinematic variables** (or a subspace thereof) such that its canonical form restricted/pushed forward equals the MHV gravity amplitude.

- Pros: directly addresses “amplitudes as canonical forms”.
- Cons: much harder; you must guess/derive the right polytope and facet inequalities.

**Recommendation:** implement **Option A** cleanly on n=4 and n=5 first (because it is a controlled theorem), then use it as a test harness for any physics-motivated map you propose.

---

## 3) Immediate engineering tasks (what to code next)

### 3.1 Make Phase G reproducible without hidden imports
`forest_polytope.py` currently imports `src.chy_oracle.forest_enumerator`. Replace that with a local enumerator or vendor the dependency, so Phase G can run standalone.

Deliverable:
- `forest_polytope.py` should run from repo root with `sage -python forest_polytope.py` (or a small wrapper), without `src/` imports.

### 3.2 Fix the lattice coordinate map in `toric.py`
`get_toric_exponents()` solves `(v - v0) = u * B_ambient` via `basis_matrix.solve_left(delta)`.
You need **guaranteed integrality**:
- Assert `u` is integral for every vertex; if not, compute a **saturated** lattice basis (Smith/Hermite normal form approach) and redo.
- Add a unit test: round-trip `v0 + u*B == v`.

Deliverables:
- `toric.py::compute_lattice_basis(..., saturate=True)`
- `tests/test_toric_roundtrip.sage` (or `.py` using Sage)

### 3.3 Implement the *actual* polytope canonical form evaluation (small n first)
Do **not** triangulate n=6 yet.
Start with:
- n=4 (dim 2 polytope; trivial)
- n=5 (dim 6; still feasible in Sage, maybe)

Implementation plan:
1. Build polytope `P = Polyhedron(vertices=...)` in its affine span.
2. Homogenize vertices to `Z_i = (1, v_i)` in `QQ^(d+1)`.
3. For a simplex `S = {i0,...,id}` define
   \[
   \Omega_S(Y) = \frac{\det(Z_{i0},...,Z_{id})}{\prod_{k=0}^d (Y\cdot Z_{ik})}.
   \]
4. Choose a triangulation of `P` and sum simplex forms.

Deliverables:
- `posgeom/canonical_form.py`:
  - `canonical_form_from_triangulation(vertices, triangulation, Y)`
  - `triangulate_polytope(vertices)` for d≤6

Validation:
- Numerically check **facet residue recursion**: residues of Ω along a facet hyperplane should match the canonical form of that facet polytope.

### 3.4 Implement the **toric pushforward check** (Option A)
Use the toric parameterization:
- choose torus coords `t_1,...,t_d`
- map into projective space via monomials `X_u = t^u` (from `get_toric_exponents`)
- implement the ABL “pushforward of dlog” check on n=4 and n=5:
  - choose a rational chart on the image
  - compute Jacobian of coordinate change
  - verify `phi^*(Ω_polytope) = const * dlog(t_1)∧...∧dlog(t_d)` numerically at random positive `t`.

Deliverables:
- `scripts/check_pushforward_n4.py`
- `scripts/check_pushforward_n5.py`

---

## 4) **Only after Option A works**: connect to physics

Phase G currently asserts “amplitude = volume / canonical form” without proving it.

### 4.1 Minimal physics-facing test (n=4)
Goal: produce a **differential form** whose pullback to kinematics matches the known 4pt MHV gravity scattering form (or at least matches the reduced Hodges determinant times a standard measure).

Procedure:
1. Use `physics_map.eval_edge_vars_from_spinors()` to get `z_ij`.
2. Decide a candidate form on `z`-space:
   - either `dlog F(z)` or `d^E z / F(z)` or a projected version.
3. Pull it back to spinor variables and compare to your oracle amplitude.

Deliverables:
- `scripts/physics_pullback_n4.sage`
- a table of random samples showing constant ratio (or not)

### 4.2 If n=4 passes, try n=5
Same pipeline as above.

If neither n=4 nor n=5 produces a clean constant-ratio comparison, **do not** scale to n=6; iterate the form ansatz.

---

## 5) Literature triage checklist (do this *before* drafting a paper)

You should explicitly read/cite at least:
- All-minors matrix-tree theorem (rooted forest determinants)
- NSVW “tree formula” for MHV graviton amplitudes
- Hodges determinant formula
- “Graphs, determinants and gravity amplitudes” (Feng et al.)
- ABL “Positive Geometries and Canonical Forms” (toric / polytope pushforward sections)
- Recent “canonical forms from adjoints” (Gaetz) for computing canonical forms without triangulating high-d polytopes
- Work relating Newton polytopes ↔ canonical forms in physics (e.g., Symanzik polytopes / recent Newton-polytope canonical form papers)

---

## 6) Exit criteria for “ready to publish”
You’re *not* ready to publish a “new physics / new positive geometry” paper until you have:

1. A clearly stated theorem-level claim that is **not** a known matrix-tree / spanning-tree identity.
2. A fully reproducible code pipeline that:
   - constructs the geometry,
   - computes its canonical form (or a well-defined pushforward),
   - and matches the physical object (Hodges/KLT/CHY) across many random kinematics with a constant ratio.
3. A novelty section that cleanly separates:
   - known identities you used,
   - your new geometric statement,
   - and your new computational verification.

