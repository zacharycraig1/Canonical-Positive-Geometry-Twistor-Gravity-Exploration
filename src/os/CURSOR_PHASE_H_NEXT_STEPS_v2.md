# CURSOR Phase H — Next Steps After Phase G (Forest Polynomial → *Actual* Canonical-Form / Pushforward Statement)

> **Context:** Phase G established a clean computational pipeline:
> - a **forest polynomial** \(F_{n,R}(z)\) (3-rooted spanning forests of \(K_n\)),
> - a **physics map** \(z_{ij}=\frac{[ij]}{\langle ij\rangle}\langle ix\rangle\langle iy\rangle\langle jx\rangle\langle jy\rangle\),
> - and numerical checks that this evaluation matches the determinantal (Hodges/Laplacian-minor) object for small \(n\).
>
> **Phase H goal:** turn this into a **precise positive-geometry / canonical-form pushforward theorem** with **correct factorization residues** — and **position it against prior literature** so we don’t “rediscover” known spanning-tree results.

---

## 0) Don’t reinvent: what Phase G already matches in the literature (read first)

Phase G’s “forest polynomial = determinant minor” identity is (mathematically) the **All-Minors Matrix–Tree Theorem** (Chaiken) applied to a Laplacian-like matrix, and (physically) closely aligned with the known “tree/forest expansions” of MHV gravity amplitudes derived from (or equivalent to) Hodges-type determinants.

**Action (required):** before writing any “new physics” claim, add a *Related Work* section that explicitly acknowledges:
- the **tree formula** / determinant ↔ tree-sum equivalence for MHV gravity,
- matrix-tree / all-minors combinatorics underlying these expansions,
- and that the **novelty target** is the **positive-geometry canonical-form / pushforward** formulation, not the existence of a forest sum.

(See “Literature & positioning” at the end for concrete papers/keywords.)

---

## 1) Immediate technical sanity checks (make the Phase G claim airtight)

### H1.1 Extend the “physics pullback” test to \(n=6\) (and optionally \(n=7\))
You already cache polytope data for \(n=6\). Now also do the *physics* identity at \(n=6\):
- Compute \(F_{6,R}(z)\) by enumerating 3-rooted forests (this is only **108** forests for \(K_6\) with 3 fixed roots — tractable).
- Evaluate under the physics map \(z_{ij}(\lambda,\tilde\lambda)\).
- Compare to the Laplacian minor/determinant construction at the same kinematics.

**Deliverable:** `src/scripts/physics_pullback_n6.sage` producing `Ratio = 1.000000` for ≥20 random seeds.

### H1.2 Gauge / reference independence stress tests
Prove numerically (then prove symbolically if feasible) that the final physical expression is independent of:
- choice of reference spinors \(x,y\),
- choice of root set \(R\) (up to a predictable overall normalization),
- choice of negative-helicity legs \(a,b\) (with corresponding \(\langle ab\rangle^8\) prefactor).

**Deliverable:** `src/scripts/gauge_invariance_sweep.sage`:
- Fix one kinematic point; sweep many random \((x,y)\) and show invariance of the *full* MHV amplitude.

### H1.3 Replace float tolerances with exact rational checks wherever possible
Many tests currently do float comparisons. When kinematics are sampled rationally, keep everything in \( \mathbb{Q}\) and verify equality exactly.

**Deliverable:** use exact `QQ` ratios, and only float for reporting.

---

## 2) Clarify the geometry: which “canonical form” are we actually claiming?

Right now, the repo computes a polytope canonical form in *dual projective space* (triangulation evaluation), and separately evaluates the forest polynomial in *edge variables* \(z_{ij}\). These are **not yet connected**.

Phase H must decide **one** of the following precise routes (and implement it):

### Route A (recommended): **Stringy canonical forms / Newton polytope → denominator = forest polynomial**
If your goal is literally “a rational form whose denominator is \(F_{n,R}(z)\)”, this matches the “stringy canonical form” / GKZ-hypergeometric style framework built from a polynomial and its Newton polytope.

**What to do:**
1. Define the Newton polytope \(P = \mathrm{Newt}(F)\).
2. Use the known construction of a canonical (or “stringy canonical”) form/function whose singular locus includes \(F(z)=0\).
3. Show that the **\(\alpha'\to 0\)** (or analogous) limit reproduces a rational function/form with denominator \(F(z)\), and that under the physics map \(z(\lambda,\tilde\lambda)\) it matches the gravity object.

**Deliverables (Route A):**
- `src/posgeom/stringy_form.py`: implement the relevant integral representation (Aomoto/GKZ style) for small \(n\) and numerically match its rational limit to \(1/F\) (or the correct \(N/F\)).
- A short note `notes/PHASE_H_ROUTE_A.md` clarifying definitions, conventions, and what exactly is being pushed forward.

### Route B: **Toric variety canonical form \(X_P\) is (tautologically) pushforward of \(d\log t\)** — but then connect it to physics
For a toric variety, the statement “pushforward of \(d\log t\) equals the canonical form on the toric positive part” is largely a **geometry identity**. The *physics* content only arrives if you:
- specify a map from **physical kinematics** into the toric variety (or its coefficients),
- and show the pullback reproduces the *physical* pole/residue structure (factorization).

**Deliverables (Route B):**
- `src/scripts/toric_pushforward_check.py`: compute Jacobians for \(t \mapsto [t^u]\) on a chosen chart and verify the form equality in coordinates (no “nonzero value” checks).
- `notes/PHASE_H_ROUTE_B.md`: explicitly state the map from kinematics → toric coordinates (or coefficients) and why its boundary behavior should match factorization.

**Important:** current `check_pushforward_n4/n5.py` only checks “nonzero evaluation”. That is not a proof of pushforward equality; Phase H must replace it with an actual differential-form check.

---

## 3) The real physics test: factorization = residues = boundaries

Even if the forest identity is correct, the “positive geometry of gravity” claim hinges on:
- **where the poles live**, and
- whether **residues on those poles** equal lower-point amplitudes (physical factorization).

### H3.1 Compute facets of the forest polytope for \(n=6\) and interpret them
The spanning tree polytope has classical facet inequalities of “subtour elimination” type; rooted-forest variants have analogous inequalities.

**What to do:**
- Compute facet inequalities for the \(n=6\) forest polytope (you claim you cache facets).
- Group facets by symmetry/orbit type.
- Attempt to map facet equations (in exponent/toric coordinates) to candidate physical limits in \(z\)-space and then to factorization channels in kinematics.

**Deliverables:**
- `src/scripts/facet_report_n6.py` producing a human-readable summary:
  - number of facets,
  - orbit decomposition under the relevant permutation subgroup,
  - explicit representative inequalities.

### H3.2 Turn facet limits into *residue computations* of the pulled-back form
Pick a boundary/facet and drive the corresponding variable(s) to the boundary:
- in toric coords \(t\), that often means \(t_i\to 0\) or \(t_i\to \infty\),
- in kinematics, this should correspond to a factorization channel (some Mandelstam \(s_S \to 0\)).

**Deliverable:** `src/scripts/residue_factorization_n6.sage` that:
1. picks a channel \(S\),
2. constructs a one-parameter kinematic family approaching \(s_S\to 0\),
3. verifies residue of the candidate form equals product of lower-point objects with correct prefactors.

### H3.3 If facet ↔ factorization mapping fails, change the geometry (don’t force it)
If the forest polytope’s facet structure does **not** match gravity factorization, that’s a clear sign the positive geometry is not “the polytope itself” but a different object (e.g., a *pullback geometry*, a blowup, a wonderful compactification chart geometry, etc.).

**Deliverable:** a decision log `notes/FACTOR_RESIDUE_DECISION_LOG.md` with:
- which facets were tested,
- what physical limit was attempted,
- whether residue matched,
- and what modification is proposed if it fails.

---

## 4) Make the “canonical form = gravity amplitude” statement precise

By the end of Phase H, you should be able to write a theorem statement **without handwaving**:

### Minimal acceptable theorem (MHV-only)
> There exists a positive geometry \((\mathcal{X}_{n,R}^{\ge 0},\Omega_{n,R})\) and a map
> \(\Phi:\mathcal{X}_{n,R}\to \mathcal{K}_n\) (kinematic data),
> such that the pushforward/pullback of \(\Omega_{n,R}\) equals the MHV gravity amplitude form, and residues on boundaries reproduce physical factorization.

To get there:
- define \(\mathcal{X}_{n,R}\) explicitly (toric variety? hypersurface complement?).
- specify the map \(\Phi\) in coordinates (not just “physics_map.py”).
- compute **at least one nontrivial factorization residue** (e.g., \(n=6\), channel \(s_{123}\to 0\)).

---

## 5) “Frontier” directions that *could* be genuinely new (choose one)

Once Phase H locks a correct canonical-form statement for MHV, the frontier moves to:

### (A) Extend beyond MHV
Try NMHV gravity, where no simple Hodges determinant exists in the same way. A positive-geometry object here would be a major step.

### (B) Loop level (even 1-loop)
Relate Newton polytopes / toric hypersurfaces / stringy canonical forms to loop integrands or leading singularities.

### (C) Replace gauge reference dependence with intrinsic geometry
Find an intrinsic formulation where \(x,y\) are eliminated geometrically (e.g., via a natural bundle construction or quotient), not “chosen then cancels”.

---

## 6) Literature & positioning checklist (do this to avoid rediscovery)

### H6.1 Mandatory web/lit scan keywords
Search and summarize (short bullets each):
- “MHV gravity tree formula”
- “Hodges determinant matrix tree theorem”
- “spanning tree / spanning forest expansion of Hodges formula”
- “stringy canonical forms Newton polytope”
- “spanning tree polytope facets Edmonds”
- “arborescence / rooted forest polytope”

### H6.2 Write `notes/RELATED_WORK.md` that explicitly says:
- **what is already known** (tree/forest expansions of MHV gravity),
- **what your novelty claim is** (positive-geometry pushforward with correct factorization residues; or a new stringy-canonical-form interpretation tied to kinematics),
- what *new consequence/prediction* comes from your geometry (e.g., a new boundary stratification / recursion / extension).

---

## 7) Concrete TODO list for Cursor (order matters)

1. **Implement** `physics_pullback_n6.sage` and pass it on ≥20 seeds.
2. **Implement** `gauge_invariance_sweep.sage` and show numerical invariance.
3. **Replace** `check_pushforward_*` with a real differential-form pushforward/pullback check (pick Route A or B, commit to it).
4. **Facet extraction** for \(n=6\) + orbit classification + write `facet_report_n6.py`.
5. **One real residue/factorization test** on \(n=6\) in a controlled kinematic limit.
6. **Write** `notes/RELATED_WORK.md` (with explicit citations) and update Phase report to clearly separate known vs novel.
7. Only after (1–6): draft “Theorem 1” and a short paper outline.

---

## Appendix: Where the repo currently stands (for quick orientation)

- Phase G summary & claims: `PHASE_G_REPORT.md`
- Forest polynomial definition: `src/posgeom/forest_polytope.py`
- Physics map \(z_{ij}\): `src/posgeom/physics_map.py`
- Canonical form eval (triangulation): `src/posgeom/canonical_polytope.py` (currently exploratory / definition-confused; Phase H must firm this up)

