# Phase E — From Spanning-Tree Geometry to the Physical 6‑pt MHV Gravity Amplitude (Cursor Agent Instructions)

**You have finished Phase D:** the Matrix–Tree (spanning tree) object is identified, its Newton polytope matches the spanning tree polytope of \(K_6\), and weights \(w_{ij}=[ij]/\langle ij\rangle\) are positive on the “positive” region you sampled.

**Phase E goal:** turn that geometric identification into an **exact algebraic bridge** from (i) the **tree/Laplacian object** to (ii) the **Hodges determinant formula / CHY integrand**, and then isolate the operation that converts the tree object into the **physical amplitude pole structure**.

This phase is intentionally “proof + invariants”: everything should be checkable with exact rational kinematics (MomentumTwistors \(\to\) spinors), and every claimed identity should be validated by random trials + scaling tests.

---

## E0. Non‑negotiable hygiene (pre-flight)

### E0.1 Use *momentum-twistor generated* kinematics for any CHY statement
Your current `phaseD4_chy_localization_mhv.py` sketch uses random spinors, which almost never satisfy momentum conservation, and CHY localization requires it.

**Fix:** for CHY-related checks, sample via `MomentumTwistor(seed=...)`, then use its `get_lambda/get_tilde_lambda`. This enforces momentum conservation through the twistor construction (in your project).

**Action**
- Create `src/chy_oracle/kinematics_samples.py` with:
  - `sample_twistor(seed)->MomentumTwistor`
  - `sample_spinors_from_twistor(seed)->(lambdas, tilde_lambdas)`
  - `assert_generic(lambdas)` checking no \(\langle ij\rangle=0\) for any required pair.

### E0.2 Eliminate floating-point where it matters
- For **identities** and **ratios**, use `QQ`/`QQbar` and exact determinants.
- For **sign** tests and quick scanning, `RR` is fine but always add an exact cross-check.

---

## E1. Prove the *correct* Matrix–Tree expansion for the Hodges reduced determinant

Right now you computed:
\[
\Sigma_{\text{Tree}}=\sum_{T\in\text{Trees}(K_6)}\prod_{(i,j)\in T} \frac{[ij]}{\langle ij\rangle}.
\]
This is the cofactor of the plain Laplacian built from weights \(w_{ij}\), **but** Hodges’ \(\Phi\) has a *reference-spinor-dependent diagonal* which corresponds to a **conjugated Laplacian** (vertex-weighted tree sum).

### E1.1 Implement Kirchhoff (determinant minor) version of the TreeSum (no enumeration)
Enumeration is fine for \(n=6\), but determinant minors will let you scale to \(n=7\) immediately and removes combinatorial overhead.

**Implement in** `src/chy_oracle/matrix_tree.py`:

- `laplacian_from_weights(w)` returning an \(n\times n\) Laplacian \(L\) with:
  - \(L_{ij}=-w_{ij}\) for \(i\neq j\)
  - \(L_{ii}=\sum_{k\neq i} w_{ik}\)
- `tree_sum_kirchhoff(w, delete=0)` returning \(\det L^{(delete)}\) (principal minor removing row/col `delete`).

**Validation**
- For \(n=6\), confirm determinant-minor equals enumeration-based sum for 10 random kinematics samples.
- Cache the list of trees only for polytope enumeration tasks; for numeric evaluation always use Kirchhoff.

### E1.2 Implement the *Hodges reference-spinor deformation* as a weighted Laplacian
Hodges diagonal uses reference spinors \(x,y\) (or your \(\chi\)-gauge choices) and effectively introduces vertex weights:
\[
C_i = \langle i x\rangle \langle i y\rangle.
\]

Define a **weighted Laplacian**:
- off-diagonal: \(\tilde L_{ij} = - w_{ij} C_i C_j\)
- diagonal: \(\tilde L_{ii} = \sum_{k\neq i} w_{ik} C_i C_k\)

Then define:
\[
\tilde \Phi := -D^{-1} \tilde L D^{-1}, \quad D=\mathrm{diag}(C_i),
\]
so that (up to the global sign convention) \(\tilde\Phi_{ij} = w_{ij}\) for \(i\neq j\) and \(\tilde\Phi_{ii} = -\sum_{k\neq i} w_{ik}\frac{C_k}{C_i}\), matching the Hodges-style diagonal.

**Action**
- Implement `hodges_weighted_laplacian(lambdas, tilde_lambdas, x, y)` producing \(\tilde L\), \(C_i\), and \(\tilde\Phi\).

### E1.3 Empirically pin the exact normalization between:
- `detprime_phi` used in your `hodges_6pt_mhv_spinor`
- the weighted Laplacian minor \(\det \tilde L^{(r)}\)
- the plain TreeSum \(\Sigma_{\text{Tree}}\)

**Deliverable script:** `src/scripts/phaseE1_tree_vs_detprime.py`

For each random twistor seed:
1. generate `(lambdas, tilde_lambdas)` from twistor
2. pick \(x,y\) (start with \(x=\lambda_0,\ y=\lambda_1\), then try random)
3. compute:
   - `detprime = detprime_phi_spinor(...)` (whatever your internal primitive is)
   - `minorL = det(tildeL.minor(r,r))`
   - `minorPhi = det(tildePhi.minor(r,r))`
4. test **exact** candidate identities:

Candidate A (conjugation relation):
\[
\det(\tilde\Phi^{(r)}) \stackrel{?}{=} (-1)^{n-1}\frac{\det(\tilde L^{(r)})}{\prod_{i\neq r} C_i^2}.
\]

Candidate B (detprime normalization):
\[
\det'(\Phi) \stackrel{?}{=} \frac{\det(\tilde\Phi^{(r)})}{(\langle x y\rangle^2\ \text{or}\ 1)}\times (\text{your deleted-row normalization}).
\]

**How to “pin” unknown factors**
- If the ratio varies, check scaling under \(x\to \alpha x\), \(y\to \beta y\).
- Any correct identity must have predictable homogeneous scaling in \(\alpha,\beta\).
- Fit the exponent of \(\langle x y\rangle\) by brute-force: try ratios multiplied by \(\langle x y\rangle^{k}\) for small integer \(k\in[-4,4]\) and see which makes ratio constant across seeds.

**Stop condition for E1:** you have a verified identity where the ratio is exactly constant across \(\ge 50\) random samples (in `QQ`), and you’ve recorded the constant + sign conventions in the report.

---

## E2. Explain (and reproduce) the *pole-structure gap* between TreeSum and the physical amplitude

Phase D observed:
- TreeSum: order \(-1\) poles at all \(\langle ij\rangle\to 0\)
- Physical amplitude: apparently different (double cyclic, suppressed/finite elsewhere in your chart)

We need to resolve this precisely **in spinor variables** and **in your chosen twistor parameterization** (these can differ).

### E2.1 Build a clean valuation engine in *spinor space*
Create `src/scripts/phaseE2_spinor_valuations.py`:

For each divisor type (adjacent and non-adjacent \(\langle ij\rangle\to 0\)):
1. start from a twistor seed (generic)
2. **force** \(\lambda_j \to \lambda_i\) at \(\epsilon=0\) by setting:
   - \(\lambda_j(\epsilon) = \lambda_i + \epsilon\,v\) with `v` chosen so \(\langle i v\rangle\neq 0\)
3. keep all other \(\lambda_k\) fixed
4. choose \(\tilde\lambda_k\) from the underlying twistor seed (do not randomize independently)
5. evaluate:
   - physical amplitude `M = hodges_6pt_mhv_spinor(...)`
   - `TreeSum_plain` via Kirchhoff on \(w_{ij}\)
   - `TreeSum_weighted` via Kirchhoff on \(w_{ij} C_i C_j\) (the one from E1)
   - any cleared numerators you care about, e.g. \(N = M\cdot D_{\text{cyclic}}\) etc.
6. fit slopes \(k(f)=d\log|f|/d\log|\langle ij\rangle|\)

**Critical**: reject trials where any other \(\langle ab\rangle\) accidentally goes to 0 (spurious mixing).

### E2.2 Decide what is “physical” vs “chart artifact”
You now have two valuation tables:
- in spinor variables: true collinear behavior
- in twistor-space chart variables: singularities induced by reconstruction denominators

**Deliverable:** `src/scripts/phaseE2_compare_spinor_vs_twistor_poles.py`
- For the *same* deformation, compute both:
  - spinor valuation (as above)
  - twistor valuation using your earlier twistor-perturbation method

**Stop condition for E2:** a written conclusion “which poles are intrinsic in spinor variables” vs “which are artifacts of the twistor parameterization.”

---

## E3. Convert Tree Geometry into the Amplitude: identify the missing operation

Phase D conclusion: TreeSum is “pre-cleared”; the physical amplitude involves a specific clearing/residue/projection.

There are three plausible mechanisms you can test **exactly**:

### Mechanism 1 — Multiply by a universal factor and take a specific minor / projection
Hypothesis:
\[
\mathcal{M}_6 \propto \frac{\langle 01\rangle^8}{(\prod \langle i,i+1\rangle^2)}\ \det(\tilde L^{(r)}) \times (\text{simple prefactor})
\]
or some close variant involving \(C_i\).

**Action**
- Use the identity from E1 and try reconstructing \(\mathcal{M}_6\) from the weighted Laplacian minor + a Parke–Taylor factor.
- In `phaseE3_reconstruct_M_from_tree.py`, scan small integer exponents \(a,b\) in:
  \[
  M_{\text{ansatz}} = \langle 01\rangle^8\ \det(\tilde L^{(r)})\ / (\prod\langle i,i+1\rangle^a)\ / (\prod_{i<j}\langle ij\rangle^b)
  \]
  and look for constant ratio against `hodges_6pt_mhv_spinor`.

### Mechanism 2 — “Cancel non-adjacent poles” via numerator vanishing (your Phase C story)
Hypothesis: the physical numerator carries systematic factors that kill some TreeSum collinear poles.

**Action**
- Use E2 valuations to **measure** the minimal power \(p_{ij}\) such that:
  \[
  \langle ij\rangle^{p_{ij}} M \quad \text{is finite as } \langle ij\rangle\to 0
  \]
- Compare this to TreeSum’s pole power (usually 1).
- This tells you exactly how many cancellations are needed and whether they are universal (same \(p\) for all non-adjacent) or helicity-dependent.

### Mechanism 3 — CHY pushforward/residue selection
Hypothesis: TreeSum is the canonical object before imposing the MHV localization + PT\(^2\) measure, and the physical amplitude is the pushforward.

**Action**
- Make the CHY side explicit in a single script `phaseE3_chy_to_tree_bridge.py`:
  - verify \((\sigma_a-\sigma_b)\propto \frac{\langle ab\rangle}{\langle a\chi\rangle\langle b\chi\rangle}\) **using momentum-twistor kinematics**
  - express the CHY integrand factor \(\prod (\sigma_{i,i+1})^{-2}\) as \(\prod\langle i,i+1\rangle^{-2}\) times \(\chi\)-dependent gauge factor
  - show where the Laplacian/tree object enters as the reduced determinant/Jacobian
- Then check numerically: CHY-localized expression equals Hodges amplitude on the same kinematics.

**Stop condition for E3:** you have a single “master identity” that maps:
**CHY localized MHV**  \(\longleftrightarrow\)  **Hodges reduced determinant**  \(\longleftrightarrow\)  **weighted Kirchhoff tree sum**, with all prefactors fixed.

---

## E4. Geometry upgrade: make the “positive polytope” statement kinematics-native

You have: spanning tree polytope in edge-variable space (generalized permutohedron).
Now you want a statement in *kinematic space*:

### E4.1 Record the exact positive region used
Your positivity scripts currently mix:
- “moment curve twistors” (true positive Grassmannian style)
- ad-hoc spinors \((1,x_i)\), \((1,y_i)\) (not necessarily momentum conserving)

**Action**
- Standardize positivity sampling:
  - Primary: moment-curve momentum twistors \(Z_i=(1,t_i,t_i^2,t_i^3)\) with ordered rational \(t_i\)
  - Derive spinors from these using your established twistor-to-spinor map
- Re-run weight positivity checks on this canonical region and confirm:
  - \(w_{ij}>0\) for all \(i<j\)
  - TreeSum and the “hat numerator” have uniform sign

Deliver `src/scripts/phaseE4_positive_region_cert.py`.

### E4.2 Identify the precise “canonical form candidate”
Once E3 fixes the prefactors, define your canonical candidate:
\[
\Omega := (\text{prefactor})\times \det(\tilde L^{(r)})\ d^{14}(\text{coords})
\]
(or whatever your derived expression is), and verify:
- log singularities only on the expected boundaries
- uniform sign on the positive region
- correct little-group scaling

---

## E5. Immediate “extend to \(n=7\)” sanity test (optional but high value)

Once you have Kirchhoff determinants, \(n=7\) becomes feasible:
- \(K_7\) has \(7^{5}=16807\) trees (enumeration is heavy but determinant is easy).

**Action**
- Add `n` support to matrix-tree routines and hodges-weighted-Laplacian routines.
- Run:
  - positivity of weights on the positive region
  - symmetry checks in the appropriate helicity-preserving subgroup
  - a small valuation table on a few boundaries
- You do **not** need Newton polytopes for \(n=7\) yet.

Deliver `src/scripts/phaseE5_n7_smoke_test.py`.

---

## Outputs required at the end of Phase E

Create/Update `src/scripts/PHASE_E_REPORT.md` with:

1. **Exact identity**: Hodges reduced determinant equals a weighted Kirchhoff tree sum (with all prefactors fixed).
2. **Pole table**: spinor-intrinsic vs chart-artifact poles (adjacent and non-adjacent).
3. **Reconstruction**: explicit formula linking TreeSum geometry to the physical amplitude (CHY/Hodges normalization).
4. **Positivity statement**: precise region definition + verified uniform sign.
5. (Optional) \(n=7\) smoke test results.

---

## Minimal “do this next” task list (agent should start here)

1. **E1.1–E1.3**: implement weighted Laplacian tree sum and match it exactly to your Hodges `detprime` primitive.
2. **E0 fix**: redo CHY localization identity using twistor-generated momentum-conserving kinematics.
3. **E2 valuations**: produce a single clean table of pole orders in spinor variables for `TreeSum_plain`, `TreeSum_weighted`, and `M`.
4. Use those to finalize **E3 master identity** (prefactors + which “clearing” produces the physical amplitude).

If any ratio fails to be constant, do not handwave—log the seed, dump all intermediate scalars (brackets, \(C_i\), minors) and add a “counterexample” section in the report.
