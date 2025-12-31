# Cursor Next Steps After Phase H (Frontier Track)

**North Star:** Don’t “rediscover Hodges/tree formulas.” Instead, either  
(1) produce a *new* positive-geometry / pushforward statement for gravity MHV, or  
(2) produce a *new obstruction/no‑go* explaining why such a geometry cannot exist in the usual “log-canonical form” sense.

This document is written so Cursor can execute immediately.

---

## 1) Reality check: what Phase G/H *really* established

### Solid progress
- You have an end-to-end pipeline that:
  1) builds the 3‑rooted forest polynomial F_{n,R}(z) on K_n,
  2) maps spinor-helicity kinematics to edge variables z_{ij} = ([ij]/<ij>) C_i C_j,
  3) compares this to a determinantal / “oracle” MHV gravity object for n=4,5,6.

- Phase H extends the check to n=6 and includes a reference-spinor sweep.

### What is **not** yet solid / publishable
- The n=6 check is described as “exact,” but the script currently uses a float ratio/tolerance. Convert it to an **exact cleared-denominator equality** (see §3).
- The “no multi-particle poles” claim is **not established**. The current script does not implement a controlled factorization limit (see §4).

---

## 2) Novelty boundary: what is likely “already known” vs what could be new

### Almost certainly already known (do NOT claim novelty here)
1) Hodges determinant formula for MHV gravity (2012).
2) Spanning-tree / rooted-forest expansions of MHV gravity amplitudes (e.g. Nguyen–Spradlin–Volovich–Wen).
3) Matrix-tree theorem / all-minors matrix-tree theorem interpretation that converts determinant minors into rooted forests.

**Cursor action:** update `notes/RELATED_WORK.md` to explicitly state:
- “We rely on Hodges + matrix-tree theorems for the forest expansion; the novelty is *not* the forest sum itself.”

### Potentially new (the only real “frontier” path)
A *precise* positive-geometry statement, such as:
- a geometry X and map φ such that **pushforward of a log form** on X equals the gravity object (amplitude or scattering form), **or**
- a canonical *function* of a polytope/toric variety whose evaluation equals the amplitude (with correct pole/soft structure),
**or**
- a proof that such a geometry cannot exist under reasonable axioms (log singularities, boundary factorization, positivity), i.e. a no‑go.

---

## 3) Immediate code fix: make n=6 verification truly exact

File: `src/scripts/physics_pullback_n6.sage`

### Replace float ratio check with exact equality
Implement:
- `lhs = M_amp * norm_factor * prod_C_sq`
- `rhs = sign * h_factor * F_val`
- Check `lhs == rhs` (Sage exact) for each trial.

Also expand coverage:
- iterate over **all** root triples R ⊂ {0..5} of size 3 (20 choices),
- randomize reference spinors (x,y) and verify invariance *exactly*,
- add a permutation test (rename labels consistently and re-check equality).

**Deliverable:** `src/tests/test_forest_identity_n4_n5_n6.sage`  
One command, returns nonzero exit status on failure.

---

## 4) Fix the physics: implement a real factorization / pole test

**You cannot claim anything about “no s_{ijk} poles”** until you do an actual limit.

Create: `src/scripts/factorization_limit_n6_s123.sage`

### Method (BCFW-style)
1) Generate generic spinor-helicity data (or momentum twistors) obeying momentum conservation.
2) Apply a BCFW shift on legs (i,j) to produce a 1-parameter family.
3) Solve for z=z* where P_{123}(z)^2 = 0 (a true 3|3 factorization channel).
4) Evaluate M_6(z) near z* and check scaling:
   - does M_6 ~ 1/(z-z*) ?
   - does the residue match the product of lower-point amplitudes?

**Outputs:**
- table of (z, P^2(z), M(z), (z-z*)M(z))
- optional plot (but keep it text-first for reproducibility)

**If factorization works:** great (aligns with physical expectations).  
**If something subtle occurs:** document carefully; that can itself be “frontier” (e.g., zeros cancelling poles).

---

## 5) Canonical form code: validate on known polytopes before trusting it on forest polytopes

File: `src/posgeom/canonical_polytope.py`

### Must-do: correctness suite
Create `src/tests/test_canonical_polytope_known_cases.py` that checks:

1) Simplex in d dimensions:
   - canonical evaluation should match the closed-form expression exactly.
2) Square/cube:
   - triangulation dependence must cancel correctly; orientation/sign matters.
3) Remove any ad-hoc absolute values or “always positive” determinant hacks.

Only after these pass should you use the routine for the forest polytope.

---

## 6) The frontier track: a *precise conjecture* and small-n pushforward computation

Create: `notes/PUSHFORWARD_CONJECTURE.md`

Write a single clean conjecture:
- Source geometry X (simplex? toric positive part of X_P? something else)
- Form Ω_X (log form)
- Map φ (explicit)
- Target object (amplitude? scattering form? canonical function?)
- Equality statement (including overall constant and domain assumptions)

### Then do n=4 explicitly
Create: `src/scripts/pushforward_n4_exact.sage`

For n=4 you should be able to compute:
- Jacobian exactly,
- pushforward exactly,
- and compare to known MHV gravity.

If n=4 fails, fix the conjecture before moving to n=5.

---

## 7) Paper positioning (so we don’t “rediscover”)
Create: `notes/PAPER_POSITIONING.md` with:

- What is known:
  - Hodges determinant, tree/forest expansions, matrix-tree theorem route.
- What is new (only one of):
  1) pushforward/canonical-object statement (with evidence),
  2) obstruction/no‑go theorem (with assumptions),
  3) new structure theorem about zeros/soft behavior encoded by polytope facets.

---

## 8) Minimal next deliverables (what Cursor should output next)

1) `src/tests/test_forest_identity_n4_n5_n6.sage`  (exact, many roots, many references)
2) `src/scripts/factorization_limit_n6_s123.sage`  (real pole test)
3) `src/tests/test_canonical_polytope_known_cases.py`  (canonical form correctness)
4) `notes/PUSHFORWARD_CONJECTURE.md`
5) Updated `notes/RELATED_WORK.md` and new `notes/PAPER_POSITIONING.md`

---

## 9) If you want a genuine “new physics” chance
The best bets are:

- **A new pushforward geometry** that naturally explains:
  - reference invariance,
  - soft/collinear limits,
  - factorization residues,
  while **not** being equivalent to an already-known determinant identity.

- Or a **no-go theorem**: gravity MHV cannot come from a positive geometry with purely log singularities + standard boundary factorization axioms. That would be legitimately valuable and publishable if done cleanly.

Either way, correctness (exact tests + real factorization limits) comes first.
