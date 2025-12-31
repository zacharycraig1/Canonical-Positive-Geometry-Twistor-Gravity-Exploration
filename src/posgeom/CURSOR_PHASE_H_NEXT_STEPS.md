# Cursor Phase H — From “Forest Polynomial = Gravity” to a **real canonical-form pushforward** (and novelty check)

> **Audience:** Cursor agent / repo maintainer  
> **Goal:** Turn Phase G’s computational identity into a *precise* positive-geometry statement (pushforward + factorization residues), while avoiding accidental rediscovery.

---

## 0) Executive reality check (what Phase G actually established vs. what it didn’t)

### Solid (and worth keeping)
- You implemented an explicit *edge-variable evaluation map*  
  \[
  z_{ij}=\frac{[ij]}{\langle ij\rangle}\,\langle i x\rangle\langle i y\rangle\,\langle j x\rangle\langle j y\rangle
  \]
  as in `physics_map.py`. fileciteturn7file10L4-L45
- You constructed the **rooted-forest generating polynomial** by brute-force enumeration in `forest_polytope.py`. fileciteturn7file9L4-L53
- You confirmed (numerically, for n=4,5) that evaluating the rooted-forest polynomial equals a reduced Laplacian minor / Hodges-type determinant (as stated in the report). fileciteturn7file1L20-L23

### Not yet established (blockers to “ready to publish”)
- Your “pushforward check” does **not** check a pushforward identity. It only evaluates *one* polytope canonical-form expression at a single random dual point and prints “non-zero”. fileciteturn7file14L32-L51
- The canonical-form evaluator is not robust enough to support claims:
  - It doesn’t enforce that each simplex has size `d+1` (the code hits a `pass` and continues). fileciteturn7file12L74-L83
  - It uses raw determinants without consistent orientation handling; even for correct triangulations this can scramble global sign/cancellations. fileciteturn7file13L48-L59
- Conceptual mismatch is openly noted inside `canonical_polytope.py`: you’re unsure whether the polytope lives in the right space for the “dual form” formula to correspond to a *physics* object. fileciteturn7file4L9-L24

**Bottom line:** Phase G is a good *infrastructure milestone* and reproduces known determinant ↔ forest-sum structure, but it is not yet a publishable new positive-geometry pushforward statement.

---

## 1) Novelty check (what is almost certainly “already known”)

The **core identity** “MHV gravity = sum over trees/forests = determinant minor” is established:

- **Tree formula for MHV gravity (forest/tree sum)**: Nguyen–Spradlin–Volovich–Wen (2009). citeturn3search1  
- **Determinant expressions / Hodges-type forms**: Hodges (2011). citeturn3search0  
- **Matrix-tree theorem link** between graph sums and determinants in gravity amplitudes: Feng & He (2012). citeturn4view2  
- **All-minors matrix-tree theorem** (forests with specified roots ↔ Laplacian minors): Chaiken (1982). citeturn3search2  

Combinatorics/geometry side:
- “Newton polytope = matroid/base polytope” for basis generating polynomials is standard; see e.g. Kummer et al. citeturn3search7  

**Implication:** A paper whose main theorem is *only* “MHV gravity equals a rooted-forest polynomial evaluated on z_ij” will read as **repackaging known results**, unless you add a genuinely new geometric pushforward/residue/factorization structure.

---

## 2) What *could* still be new (and worth pursuing)

Plausibly new if you can supply **both**:
1. A **positive geometry** with a clearly defined map such that its **canonical form** matches the gravity object after pushforward/pullback; and  
2. Boundaries/residues that match **physical factorization channels**.

There is precedent: Lam explains the associahedron story as an **algebraic moment map** of a toric variety whose pushforward gives a polytope canonical form, citing the AHBL pushforward theorem. citeturn4view0turn4view1

---

## 3) Phase H deliverables (what Cursor should build next)

### H1. Promote the identity into a clean “known background” theorem
Create `docs/theorem_inventory.md`:
- Precise statement actually verified numerically (and how it matches Laplacian minors).
- Proof sketch via **Chaiken all-minors matrix-tree theorem**. citeturn3search2  
- Related-work section: NSVW, Hodges, Feng–He. citeturn3search1turn3search0turn4view2  

### H2. Make the canonical-form code trustworthy
Modify `src/posgeom/canonical_polytope.py`:
1. **Triangulation sanity:** assert simplex size = `d+1`; if not, raise (no silent `pass`). fileciteturn7file12L74-L83  
2. **Orientation convention:** either use `abs(det)` + global sign, or orient simplices consistently. fileciteturn7file13L48-L59  
3. **Triangulation invariance test:** evaluate Omega(W) for two triangulations at ~20 random W and confirm equality.

### H3. Implement the *actual* toric pushforward (A.6 for real)
Add `src/posgeom/pushforward_toric.py` and define the full diagram (spaces + maps).

#### Option H3A (AHBL/Lam-style): algebraic moment map → polytope form
1. Build projective toric variety \(X_P\) from your polytope \(P\) via the monomial map  
   \[
   \phi:(\mathbb{C}^*)^d\to \mathbb{P}^{V-1},\quad t\mapsto [t^{u_1}:\cdots:t^{u_V}].
   \]
2. Implement the **algebraic moment map** \(\mu\) used in AHBL/Lam-style pushforwards. citeturn4view0turn4view1  
3. Numerically verify pushforward: compare \((\mu\circ\phi)_*(d\log t_1\wedge\cdots\wedge d\log t_d)\) to \,\(\Omega(P)\).

**Acceptance test:** n=4 and n=5 match (up to global sign) on many random positive samples.

> If you cannot define \(\mu\) concretely, you do not yet have a pushforward statement.

#### Option H3B (more physics-native): build the polytope in kinematic space
Route B: construct a polytope whose facets are physical kinematic hyperplanes; show its canonical form pulls back to gravity (à la associahedron scattering forms). citeturn4view0

### H4. Residue / factorization matching (the “physics test”)
After H3:
- Enumerate facets/toric divisors.
- For each, compute residues and match to products of lower-point objects.
Lam’s method: prove equality by matching poles + residues. citeturn4view0turn4view1

### H5. Scale to n=6 without pain
Replace spanning-tree enumeration with a direct rooted-forest generator for \(K_n\) (Prüfer-like). Your current approach will not scale. fileciteturn7file9L16-L53  

**Acceptance test:** for n=6, roots R=[0,1,2], you get exactly 108 forests in <1s.

---

## 4) Publish/no-publish gating criteria

**Not ready to publish if:**
- No explicit map and no verified pushforward equality.
- No residue/factorization checks.
- No explicit related-work section acknowledging NSVW/Hodges/Feng–He/Chaiken.

**Potentially publishable if:**
- You prove a new theorem of the form  
  “There exists a positive geometry \(\mathcal{G}_n\) with map \(f\) such that \(f_*\Omega(\mathcal{G}_n)=\Omega_{\text{grav},n}\)”  
  and verify poles+residues at n=6 (and ideally general n).

---

## 5) Concrete next files Cursor should add / modify

### Create
- `docs/theorem_inventory.md`
- `src/posgeom/pushforward_toric.py`
- `src/tests/test_canonical_form_polytope.py`
- `src/tests/test_pushforward_n4_n5.py`
- `src/scripts/physics_pullback_n6.sage` (extend identity check to n=6)

### Modify
- `src/posgeom/canonical_polytope.py` (no silent passes; orientation; invariance tests)
- `src/posgeom/forest_polytope.py` (direct forest generator)

---

## 6) Repo bibliography (keep these citations)
- Nguyen–Spradlin–Volovich–Wen (2009). citeturn3search1  
- Hodges (2011). citeturn3search0  
- Feng & He (2012). citeturn4view2  
- Chaiken (1982). citeturn3search2  
- Arkani-Hamed–Bai–Lam (2017). citeturn4view1  
- Lam, “Moduli spaces in positive geometry” (moment maps + pushforward). citeturn4view0  
