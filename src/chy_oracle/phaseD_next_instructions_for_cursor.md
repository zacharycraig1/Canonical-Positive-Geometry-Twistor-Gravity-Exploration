# Phase D — From “Hat Numerator Positivity” to an Explicit Positive Geometry (CHY-native Northstar)

You’ve completed Phase C3: a fast, exact **spinor-based Hodges/CHY oracle**, boundary **valuation tooling**, rigorous **helicity-preserving symmetry checks**, and a **positivity signal** (uniform sign of \(\hat N\) on the positive region via moment-curve kinematics).

**Phase D goal:** convert that positivity signal into an **explicit, geometry-native object**: a *named* positive geometry (polytope / matroid polytope / cluster polytope / wonderful compactification chart object), plus a **proof-grade bridge** from CHY/Hodges to that geometry.

This doc gives Cursor a concrete roadmap with implementable milestones and “what success looks like.”

---

## D0 — Unify conventions (one “master cleared object”)

You currently have (at least) two “clearing denominators” in play:

- **Twistor-space discovery (Phase B):**  
  \[
  D_{\text{found}}(Z)=\frac{\Big(\prod_{i<j}\langle ij\rangle\Big)\Big(\prod_{\text{cyc}}\langle i,i{+}1\rangle\Big)}{\langle 01\rangle^2}
  \]
  and \(N_{\text{found}} = M_{\text{CHY}} \cdot D_{\text{found}}\) was polynomial on lines \(Z(t)\).

- **Spinor-space clearing used in C3:**  
  \[
  D_{\text{cyc}}(\lambda)=\Big(\prod_{\text{cyc}}\langle i,i{+}1\rangle\Big)^2,\qquad
  N_{\text{cyc}} := M \cdot D_{\text{cyc}},\qquad
  \hat N := \frac{N_{\text{cyc}}}{\langle 01\rangle^8}.
  \]

**Task D0.1 (required):** write a short “conventions ledger” in `src/scripts/PHASE_D_CONVENTIONS.md` that pins down:

1. What exactly `M` is in each pipeline (CHY integrand evaluation vs Hodges determinant formula vs any extra Jacobians from parameterizations).
2. A precise statement of what space each denominator lives in:
   - momentum twistors \(Z\) (Gr(4,6) / momentum-twistor chart),
   - spinors \((\lambda,\tilde\lambda)\) (on-shell phase space modulo little group),
   - kinematic invariants \(s_{ij}\) (kinematic space).
3. The mapping relationships: which extra factors are *pure Jacobians / chart artifacts* vs true physical poles.

**Success =** you can state (and verify numerically on random samples) an identity of the form
\[
M\cdot D_{\text{found}}(Z)\;\propto\;(M\cdot D_{\text{cyc}}(\lambda))\cdot J(\lambda,\tilde\lambda)
\]
for some explicitly computable “map Jacobian” \(J\), or at least isolate the mismatch to a known chart artifact.

> Don’t over-engineer D0: we just need a stable dictionary so Phase D claims are not “representation-dependent.”

---

## D1 — Exploit the Laplacian / Matrix–Tree theorem (this is the best geometry-native lever)

The Hodges reduced determinant is a **Laplacian minor**, hence equal to a **Kirchhoff (spanning tree) polynomial** by the Matrix–Tree theorem. This is the most “positive geometry adjacent” structure you have access to *right now*.

### D1.1 Implement a matrix-tree expansion oracle
Create `src/chy_oracle/matrix_tree.py` with:

- `hodges_minor_matrix_tree(lambdas, tilde_lambdas, delete=(0,1,2), diag_gauge="hodges")`
  - builds weights \(w_{ij}=[ij]/\langle ij\rangle\) (for \(i\ne j\)),
  - returns \(\det'(\Phi)\) via **sum over spanning trees** on \(K_6\) (or \(K_n\) later),
  - includes caching of spanning trees of \(K_6\) (there are \(6^{6-2}=1296\) trees; trivial to precompute once).

- `compare_det_vs_tree(ntrials=50)`
  - confirms exact equality with `det`-based reduced determinant for random seeds.

**Success =** tree-sum equals determinant oracle exactly over \(QQ\), 50/50 trials.

### D1.2 Factor out denominators to get a *true polynomial with positive coefficients*
A tree term is \(\prod_{(i,j)\in T} \frac{[ij]}{\langle ij\rangle}\).  
This is not polynomial as-is. But your observed “uniform sign on the positive region” strongly suggests a *canonical clearing* exists that converts the tree-sum into a **polynomial with all coefficients \(+1\)** (or fixed sign) in **some variables**.

Implement in `src/scripts/phaseD1_tree_clearing.py`:

1. Choose a candidate global factor \(F(\lambda)\) (start with your empirically good \(D_{\text{cyc}}/\langle 01\rangle^8\), but also test other natural candidates like \(\prod_{i<j}\langle ij\rangle\)).
2. Define
   \[
   P_T := F(\lambda)\cdot \det'(\Phi)
   \]
   and attempt to confirm polynomiality by:
   - line sampling in \(\lambda(t),\tilde\lambda(t)\) using your existing rational reconstruction harness,
   - or direct “pole cancellation” checks: verify no remaining negative valuation in \(\langle ij\rangle\) directions.

3. Once \(P_T\) is polynomial-like, **extract monomial supports**:
   - Evaluate \(P_T\) on a family where each bracket is assigned a symbolic monomial weight (see D2.2), or
   - Use repeated random specialization and sparse reconstruction (finite fields can help here).

**Success =** you produce a *closed-form statement*:
- either “\(P_T\) is polynomial in \(\langle ij\rangle\) with integer coefficients,” or
- (better) “\(P_T\) is the evaluation of a **Kirchhoff polynomial** (tree polynomial) in positive edge variables.”

---

## D2 — Turn “uniform sign” into a theorem (sign patterns + tree structure)

Your C3 result: \(\hat N\) is uniformly negative on moment-curve positive kinematics.

Tree-sum interpretation gives a route to *prove* this:
- If every tree term has the same sign on the positive region, the sum is sign-definite.

### D2.1 Build a sign table for weights \(w_{ij}=[ij]/\langle ij\rangle\)
Write `src/scripts/phaseD2_weight_signs.py`:

- generate 200–1000 positive moment-curve samples,
- compute the sign of each:
  - \(\langle ij\rangle\)
  - \([ij]\)
  - \(w_{ij}\)
- summarize:
  - whether signs of \(w_{ij}\) are uniform,
  - or whether a simple parity rule exists (e.g. sign depends only on \(i+j\) mod 2).

**Success =** you discover a deterministic sign rule for \(w_{ij}\) on the positive region.

### D2.2 Prove sign-definiteness of the tree sum from the rule
If the sign rule implies \(\prod_{(i,j)\in T} w_{ij}\) has the same sign for every tree \(T\), then:
\(\det'(\Phi)\) has fixed sign; multiplying by a manifestly positive bracket factor yields \(\hat N\) fixed sign.

Write `src/scripts/phaseD2_tree_term_signs.py`:

- For each sample, compute signs of all 1296 tree terms (fast) and verify “all equal.”
- Store counterexamples if any.

**Success =** 0 counterexamples over 1000 samples; then write a proof sketch in `PHASE_D_REPORT.md`.

---

## D3 — Identify the polytope: Newton polytope = spanning tree polytope (matroid base polytope)

If you succeed in D1–D2, you’re likely sitting on a known combinatorial geometry:

- The Kirchhoff polynomial of a graph has Newton polytope equal to the **spanning tree polytope** (a **graphic matroid base polytope**).
- For the complete graph \(K_6\), this polytope is a classic **generalized permutohedron** object.

### D3.1 Compute Newton polytope support (exponent vectors)
In `src/scripts/phaseD3_newton_polytope.py`:

- represent each monomial by an exponent vector over “edge variables” \(x_{ij}\) (one per edge),
- for each tree term, exponent vector is 1 on its edges, 0 elsewhere,
- build the convex hull (Sage `Polyhedron`) and compute:
  - dimension,
  - facet inequalities,
  - vertex count (should match #trees for extreme points, but check).

**Success =** a recognized facet structure consistent with spanning tree polytopes (cut constraints).

### D3.2 Compare to known descriptions
Write short notes in `PHASE_D_REPORT.md` explaining:
- how your \(\hat N\) maps to the Kirchhoff polynomial,
- how the Newton polytope matches the spanning tree polytope,
- why “associahedron + permutohedron” intuition appeared earlier.

This is the point where the project becomes “geometry-native” in a way a math/physics audience recognizes immediately.

---

## D4 — CHY-native geometric statement on \(M_{0,6}\) (worldsheet geometry)

You now have three spaces in play: \(M_{0,6}\) (punctures), kinematic space, and twistor/spinor space.

### D4.1 Localize the CHY formula directly to the MHV solution
Create `src/scripts/phaseD4_chy_localization_mhv.py` that:

- takes random rational kinematics,
- builds the 4D MHV solution \(\sigma_a\) you already implement (`sigma_mhv`),
- evaluates the CHY Jacobian and integrand pieces,
- outputs a factorized expression showing where:
  - Parke–Taylor (YM) and Hodges determinant (gravity) arise,
  - and what factors correspond to worldsheet boundary divisors \(\sigma_i=\sigma_j\).

**Success =** a clear factorization log that aligns “worldsheet divisor \(\leftrightarrow\) bracket pole” and explains why certain poles are doubled in some parameterizations.

---

## D5 — Scale up a tiny bit: n=7 smoke test (pattern check, not a full rebuild)

Do a “small generalization” purely as a robustness check:

- implement `hodges_7pt_mhv_spinor` (same approach),
- repeat:
  - symmetry checks (helicity-preserving subgroup),
  - positivity probe on moment curve,
  - quick weight-sign table.

**Success =** the sign-definiteness persists, suggesting a genuine positive-geometry pattern rather than a 6-pt accident.

---

## Deliverables for Phase D (what to commit)

1. `src/scripts/PHASE_D_CONVENTIONS.md` — the dictionary between the pipelines.
2. `src/chy_oracle/matrix_tree.py` — tree-sum oracle + cache.
3. `src/scripts/phaseD1_tree_clearing.py` — shows how \(\hat N\) emerges as a polynomial / evaluated Kirchhoff polynomial.
4. `src/scripts/phaseD2_weight_signs.py` and `phaseD2_tree_term_signs.py` — positivity theorem machinery.
5. `src/scripts/phaseD3_newton_polytope.py` — spanning tree polytope computation.
6. `src/scripts/phaseD4_chy_localization_mhv.py` — explicit CHY localization narrative.
7. `src/scripts/PHASE_D_REPORT.md` — readable “paper-style” report tying all of the above together.

---

## Practical engineering notes (keep runs fast and stable)

- **Precompute spanning trees once** (1296 for \(K_6\)); store as list of edge lists in `.sobj` or `.json`.
- **Avoid repeated `load(hodges.sage)`** inside scripts; centralize imports.
- For sign/positivity probes, stay in **exact rationals** as long as possible; if you need speed, switch to **mod prime finite field** sampling for structure discovery, then lift via reconstruction.
- In valuation scripts, reject seeds where any neighbor bracket is exactly 0 before perturbation (cheap guard).

---

## “What’s the breakthrough claim?” (keep it honest and sharp)

If D1–D3 succeed, you can claim:

> The 6‑pt MHV gravity Hodges determinant (after a natural clearing) evaluates a **Kirchhoff / spanning-tree polynomial** with **sign-definite weights** on the positive kinematic region, and its Newton polytope is the **spanning tree polytope of \(K_6\)** (a generalized permutohedron). This provides a concrete, geometry-native “positive” object underlying the amplitude, directly tied to CHY localization.

That’s a crisp, defensible geometry-native statement — and it naturally generalizes.

---

## Immediate next command for Cursor

1. Implement `matrix_tree.py` + cache spanning trees for \(K_6\).
2. Run equality check vs reduced determinant oracle (50 trials).
3. Run `phaseD2_weight_signs.py` on 500 moment-curve samples.
4. If sign rules look clean, run `phaseD2_tree_term_signs.py` to confirm “all trees same sign.”
