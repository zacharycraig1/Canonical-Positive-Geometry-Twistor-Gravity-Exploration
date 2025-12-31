# Phase C2 — Robust Channel Factorization (BCFW) — Fixes + Optimization + Northstar

**Context:** You are testing the 6‑pt MHV gravity amplitude in a CHY/Hodges oracle framework.  
**Goal:** Verify (and correctly interpret) behavior as multi‑particle channels (e.g. \(s_{123}\)) approach zero, *without* getting tricked by parameterization poles, and make the scripts fast + stable for large automated runs.

---

## 1) What the latest results mean (physics interpretation)

Your BCFW scan is consistently finding exact \(z_\*\) values where \(s_{123}(z_\*) = 0\), and the amplitude stays **finite** with \(\alpha \approx 0\).

That is **plausible and expected** for **6‑pt MHV gravity with only two negative helicities**:

- In the factorization limit,
  \[
  \mathcal{M}_6 \xrightarrow{s_{123}\to 0} \sum_h \mathcal{M}_L(1,2,3,P^h)\,\frac{1}{s_{123}}\,\mathcal{M}_R(-P^{-h},4,5,6).
  \]
- Each side is a **4‑pt** gravity amplitude.
- The right side has external \((+,+,+)\) and the internal line provides at most one \((- )\) → **single-minus or all-plus tree gravity amplitudes vanish**, so the residue is zero.
- Therefore there is **no \(1/s_{123}\)** singularity: \(\mathcal{M}_6\sim s^0\) is expected.

✅ **Conclusion:** “No channel pole” for \(s_{123}\) is *not* necessarily an error; it’s a helicity selection rule story.

---

## 2) Fix the Sage `roots()` keyword error

Your log shows:

> `roots() got unexpected keyword argument(s): dict_keys(['multiplicity'])`

In Sage, the keyword is **`multiplicities`**, not `multiplicity`.

Use one of:

```python
roots = P.roots(ring=QQbar, multiplicities=False)  # returns a list of roots only
# OR
roots = P.roots(ring=QQbar)                        # default returns (root, mult) pairs
```

If you want multiplicities explicitly:

```python
roots = P.roots(ring=QQbar, multiplicities=True)   # (root, mult)
```

---

## 3) Docker: the correct way to run Sage scripts for this project

Two common confusions:
- `sage` “CommandNotFoundException” in PowerShell means you ran outside the container/WSL.
- A container “exits immediately” can be **normal** if it ran your command and finished (exit code 0).

### Run the script (PowerShell, from the project folder)

```powershell
docker run --rm -it -v ${PWD}:/workspace -w /workspace sage-cursor:latest `
  sage -python src/scripts/phaseC2_factorization_robust.py
```

### Open an interactive shell inside the container

```powershell
docker run --rm -it -v ${PWD}:/workspace -w /workspace sage-cursor:latest bash
```

Then you can run:

```bash
sage -python src/scripts/phaseC2_factorization_robust.py
```

---

## 4) Script hygiene fixes (avoid Cursor merge corruption)

### A) Fix broken imports (you had `import sysfrom sage.all import *import os`)
Replace the import block with:

```python
import sys, os, json, time
from sage.all import *
```

### B) Fix `NameError: get_s123 is not defined`
Ensure helper functions are defined **above** any function that calls them, especially above:

- `run_robust_factorization_check_bcfw()`

Order should be:

1. imports
2. `sys.path.append(...)`
3. `load(hodges.sage)` (or imports)
4. helpers: `ang_bracket`, `get_s123`, `bcfw_shift_spinors`, `solve_s123_bcfw`, etc.
5. main run function
6. `if __name__ == "__main__": ...`

---

## 5) Biggest performance + stability optimization (recommended)

Right now you do:

**BCFW shift (spinors)** → **reconstruct momentum twistors** via `spinors_to_twistors()` → evaluate `hodges_6pt_mhv(twistor)`.

This is expensive and sometimes singular (e.g. `domain_violation_angle_bracket_offdiag`).

### ✅ Replace with “spinor-only Hodges” for Phase C2 scaling tests

Implement:

- `hodges_6pt_mhv_spinor(lambdas, tilde_lambdas, neg=(0,1))`

Core idea:
- build Hodges matrix \(\Phi\) with off-diagonal
  \[
  \Phi_{ij} = \frac{[ij]}{\langle ij\rangle},\quad i\ne j
  \]
- set diagonal by the Hodges prescription (sum rule with reference spinors / momentum conservation form)
- take reduced determinant
- multiply by the helicity prefactor \(\langle 01\rangle^8\) (or your chosen convention)

Then in the scaling loop:
- compute \(s_{123}\) directly from spinors
- compute \(\mathcal{M}\) directly from the spinor Hodges implementation
- **do not** map back to momentum twistors for this test

Expected benefits:
- Much faster per trial
- Removes most domain violations
- Makes the “\(\alpha\)” scaling estimates more trustworthy

Use momentum twistors only when you specifically need twistor-geometry statements about \(D(Z)\), \(N(Z)\), etc.

---

## 6) Northstar checklist for the bot (next concrete milestones)

### (1) Channel classification by helicity selection rules
For each multi-particle channel (e.g. \(s_{123}, s_{234}, s_{345}, \ldots\)):

- use a BCFW shift that **actually moves** that channel
- solve exactly for \(z_\*\) such that \(s_S(z_\*)=0\)
- measure scaling:
  - does \(M\sim 1/s\)? (nonzero residue)
  - or \(M\sim s^0\)? (vanishing residue / helicity-forbidden)
- store results in JSON: channel, seed, z*, alpha

### (2) Confirm twistor singularities are only bracket boundaries
Using your discovered clearing denominator \(D(Z)\), verify on multiple random lines:
- \(N(Z)=M\cdot D(Z)\) is polynomial
- pole orders match:
  - cyclic \(\langle i,i+1\rangle\): double
  - non-cyclic \(\langle i,j\rangle\): simple
- \(\langle 01\rangle\) cancellation order matches your measured \(\langle 01\rangle^6\) factor

### (3) Geometry / positivity probe (if desired)
On positive kinematics (moment curve samples):
- evaluate sign/structure of \(\hat N = N/\langle 01\rangle^6\)
- locate the zero locus \(\hat N=0\)
- check whether it behaves like a canonical-form numerator (sign patterns, boundary behavior)

---

## 7) Quick triage if you hit issues again

- If `sage` not found → you are not running inside Sage (container/WSL). Use the docker commands above.
- If factorization search “blows up” → you’re hitting a **pole** of the parameterization; reject that interval/seed.
- If `spinors_to_twistors` fails → skip that trial or (better) switch to the spinor-only Hodges implementation.

---

**Deliverables to implement next (recommended):**
1. Patch `roots()` keyword usage (`multiplicities`).
2. Clean import block + function ordering in `phaseC2_factorization_robust.py`.
3. Add `hodges_6pt_mhv_spinor(...)` and update PhaseC2 to avoid twistor reconstruction.
4. Run channel classification for all partitions \(S\subset\{1..6\}\) with |S|=3 (unique channels at n=6).
