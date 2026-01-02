# Cursor Phase O+: From “Forest Identity” to a Real Canonical-Form Pushforward (Frontier-Focused)

**Context:** You’ve now got a *very solid* internal oracle (Hodges reduced determinant) with enforced conventions + nontrivial checks passing (deletion invariance, soft scaling, collinear scaling). fileciteturn43file1L7-L21  
**Blocker:** the independent KLT oracle still disagrees (ratio not constant), so we *cannot* claim an independent verification yet. fileciteturn43file1L16-L21  
**Bigger issue:** the “canonical-form/pushforward” statement is still not mathematically pinned down (the doc has “?” placeholders). fileciteturn43file5L45-L55

This doc gives the next steps that (a) **avoid rediscovering known results** and (b) **target the smallest publishable frontier advance**: a precise positive-geometry statement with a correct pushforward map and residue/factorization interpretation.

---

## 0) Novelty reality-check (do not overclaim)

Before you write any paper claiming “new theorem,” internalize:

### What is *already known*
1. **Hodges’ determinant formula** for MHV gravity is standard.
2. **Matrix-tree / spanning tree/forest expansions** of the Hodges determinant are known; there are explicit “tree formula” representations of MHV gravity amplitudes in the literature (e.g. NSVW tree formula; and related determinant/graph expansions by Feng–He and others).  
   **Meaning:** “MHV gravity = (some determinant) = sum over trees/forests” is *not* new.

### What could still be new (if you nail it)
- A **precise positive-geometry** whose **canonical form** (or **stringy canonical form** in the α′→0 sense) **pushes forward** to the gravity MHV amplitude under an explicit map from a positive domain (torus / moduli / kinematic region).  
- A **clean boundary/facet ↔ physical limit** dictionary for that geometry that is *not* just “the determinant expands into trees.”

Your goal is to turn your current “forest polynomial identity” results into one of the above. Anything less is “repackaging known structures.”

---

## 1) Phase O: Fix the independent oracle mismatch (KLT vs Hodges)

**Deliverable:** `RESULTS/klt_vs_hodges_certificate.json` showing constant ratio for n=4,5,6 (or a pinpointed convention bug if not).

### O1. Implement *standard* KLT with a reference formula (do not trust the current kernel)
Your current end-to-end harness reports a non-constant ratio. fileciteturn43file1L16-L21  
Assume KLT is wrong until proven otherwise.

#### O1.a Standard KLT formula (field theory)
Pick labels **1..n** (convert your 0..n-1 indexing carefully).

A common field-theory KLT form is:
\[
M_n \;=\; (-1)^{n+1}\!\!\sum_{\alpha,\beta\in S_{n-3}}
A_n(1,\alpha,n-1,n)\; S[\alpha|\beta]_{1}\; A_n(1,\beta,n,n-1),
\]
where \(A_n\) are color-ordered Yang–Mills tree amplitudes and \(S[\alpha|\beta]_1\) is the momentum kernel.

**Momentum kernel (algorithmic definition):**  
Let \(\alpha = (\alpha_1,\ldots,\alpha_{n-3})\), \(\beta=(\beta_1,\ldots,\beta_{n-3})\) permutations of \(\{2,\ldots,n-2\}\). Then
\[
S[\alpha|\beta]_1 \;=\;\prod_{i=1}^{n-3}
\left( s_{1\alpha_i} + \sum_{j<i} \Theta(\alpha_j,\alpha_i;\beta)\, s_{\alpha_j\alpha_i}\right),
\]
where \(\Theta(\alpha_j,\alpha_i;\beta)=1\) if \(\alpha_j\) appears **before** \(\alpha_i\) in the list \(\beta\), and \(0\) otherwise.

**Critical:** ensure your Mandelstam convention matches what KLT expects. You standardized
\(s_{ij} = \langle ij\rangle[ji]\). fileciteturn43file1L7-L10  
That’s fine, but the sign matters inside the kernel; test small n first.

#### O1.b Parke–Taylor for YM MHV (your ordering must match kernel ordering)
For a cyclic ordering \((\sigma_1,\ldots,\sigma_n)\) with negative helicities on legs \(a,b\):
\[
A_n^{\text{MHV}}(\sigma) \;=\; \frac{\langle a b\rangle^4}{\langle \sigma_1\sigma_2\rangle\langle \sigma_2\sigma_3\rangle\cdots\langle \sigma_n\sigma_1\rangle}.
\]
(Use the ordering \(\sigma\) in the denominator exactly.)

#### O1.c Ground-truth tests (must pass before n=6)
1. **n=4**: Implement
\[
M_4 = -\, s_{12}\, A(1,2,3,4)\, A(1,2,4,3).
\]
2. **n=5**: Implement the 2×2 KLT sum over permutations of \(\{2,3\}\) (or \(\{2,3\}\) depending on fixed legs).
3. Only then move to **n=6**.

### O2. Add a 3rd oracle (avoid single-point failure modes)
If KLT is painful, add an independent check:
- **BCFW recursion** for gravity MHV (more work but very clean), or
- **CHY gravity** (Pfaffian-squared) evaluated numerically (harder but canonical).

**Acceptance criterion:** Hodges matches at least one independent method at random kinematics with a constant ratio (±1, or an explainable convention constant).

---

## 2) Phase P: Make the Forest Polytope definition “hard” (inequalities + facets)

Right now “Forest polytope = Newton polytope of forest polynomial” is OK, but you need an explicit *facet/inequality* model to talk about boundaries and residues.

### P1. Define the rooted-forest polytope as a face of a spanning-tree polytope
A clean way:

1. Start with graph \(G = K_n\) on vertices \([n]\).
2. Add a **super-root** vertex \(\rho\) and edges \((\rho,r)\) for each \(r\in R\).
3. Consider spanning trees of \(G' = K_n \cup \{(\rho,r)\}\) **that include all \((\rho,r)\)** edges.
4. Removing \(\rho\) and its incident edges gives a spanning forest of \(K_n\) with components rooted at \(R\).

This identifies your forest set with a **face** of the graphic matroid base polytope of \(G'\).

### P2. Inequalities (practical implementation target)
For the spanning tree polytope \(P_{\text{tree}}(G')\) in variables \(x_e\):
- \(x_e \ge 0\)
- \(\sum_{e\in E(G')} x_e = |V(G')|-1\)
- For all nonempty \(S \subset V(G')\): \(\sum_{e\in E(S)} x_e \le |S|-1\)

Then intersect with the face constraints \(x_{(\rho,r)}=1\) for \(r\in R\), and project away the \((\rho,r)\) coordinates.

**Tasks:**
- Implement `src/posgeom/forest_polytope_inequalities.py` generating these inequalities for given (n,R).
- For n=6, verify Sage facet enumeration matches (up to redundancy) your inequality list.

**Deliverable:** `RESULTS/facet_certificate_n6.json` mapping each computed facet to an inequality family.

---

## 3) Phase Q: Replace “canonical form = polynomial” with a correct canonical-form object

### Q1. Fix the canonical form evaluator (it’s currently fragile)
Your `eval_canonical_form_dual` sums signed simplex determinants from a triangulation: fileciteturn43file8L37-L63  
and uses `vol = det(...)` without an orientation guarantee. fileciteturn43file16L45-L60  
This can silently give the wrong rational function.

**Do this:**
- Add an oriented triangulation (or take absolute determinants with a consistent sign convention derived from a reference simplex).
- Add a pole-avoidance sampler for W (reject if any \(W\cdot Z_v=0\)).

**Unit tests:**
- For an actual simplex, match the closed form:
\[
\Omega_{\Delta}(W)=\frac{\det(Z_0,\dots,Z_d)}{\prod_{i=0}^{d}(W\cdot Z_i)}.
\]
- For a square (2D), compare two triangulations and verify equality.

### Q2. Decide which “canonical form” you mean (and document it)
There are *three* different things floating around; pick one and be consistent:

1. **Polytope canonical form** on dual projective space (rational function in W).
2. **Stringy canonical form** \(I_P(s;\alpha')\) whose \(\alpha'\to 0\) limit yields the polytope canonical form in an \(s\)-space.
3. **Log form on a hypersurface complement** like \(d\log F\wedge\cdots\) (not a polytope canonical form).

Your current pushforward doc literally has “?” here: fileciteturn43file5L45-L55  
so you must resolve it.

---

## 4) Phase R: The actual frontier push—prove a pushforward theorem, not an identity

Once O, P, Q are solid, attempt *one* tight publishable statement:

### R1. Candidate statement (template)
Define a positive domain \(\mathcal{D}\) (e.g., a positive torus chart / a moduli-space cell / a kinematic positive region) with canonical log form \(\Omega_{\mathcal{D}}\).

Define a map \(\Phi:\mathcal{D}\to \mathbb{P}^d\) (or kinematic space) such that:
\[
\Phi_*(\Omega_{\mathcal{D}}) \;=\; \Omega_{P_F}(\cdot)
\]
and then show that after composing with the *physics map* \( \Psi(\text{spinors})\) you get:
\[
M_n^{\text{MHV}}(\lambda,\tilde\lambda)\;=\; \left(\Psi^*\circ \Phi_*\right)(\Omega_{\mathcal{D}}).
\]

**Important:** this is where novelty can live.

### R2. Boundary dictionary (the hard physics content)
Show that each relevant physical limit corresponds to a facet/face degeneration:
- Soft limit = contraction/deletion specialization of forest polynomial + facet behavior
- Collinear limit = specific face restriction that reproduces splitting

You already have soft/collinear scaling checks passing at the amplitude level. fileciteturn43file1L10-L15  
Now you must tie them to **specific facets** from Phase P.

---

## 5) Minimal run commands (make Cursor deterministic)

Add to `AUTO_RUN.md` (or create `RUN_PHASE_O.md`):

1. **KLT sanity:**
   - `sage -python src/chy_oracle/test_klt_n4.py`
   - `sage -python src/chy_oracle/test_klt_n5.py`
   - `sage src/pushforward/end_to_end_n6.sage`

2. **Facet certification:**
   - `sage -python src/posgeom/forest_polytope_inequalities.py --n 6 --roots 0,1,2 --emit RESULTS/facets_n6.json`
   - `sage -python src/posgeom/compare_facets.py`

3. **Canonical-form unit tests:**
   - `pytest -q src/tests/test_canonical_simplex.py`
   - `pytest -q src/tests/test_canonical_square.py`

---

## 6) “Stop conditions” (avoid infinite loops)

Stop and write up if you achieve **any** of:
- KLT = Hodges (constant ratio) for n=4,5,6 on ≥50 random kinematics points.
- Facet inequalities proven + matched to computed facets for n=6.
- A coherent, non-handwavy pushforward map \(\Phi\) with at least one nontrivial Jacobian/pushforward check in a toy model (n=4,5) that survives triangulation-independence tests.

If you don’t have these, you’re not ready to claim a new positive-geometry theorem.

---

## Appendix: what *is* correct right now (and what is not)

### Correct / strong
- Mandelstam convention and Hodges deletion invariance were explicitly tested and passed. fileciteturn43file1L7-L10
- The KLT-vs-Hodges mismatch is explicitly acknowledged (good). fileciteturn43file3L23-L30

### Not yet correct/finished (must fix)
- The pushforward/canonical-form definition is still ambiguous (“?” in the doc). fileciteturn43file5L45-L55
- The canonical-form evaluation needs orientation/triangulation-independence hardening. fileciteturn43file16L45-L60
