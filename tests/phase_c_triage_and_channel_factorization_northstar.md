# Phase C Triage + North Star Update (Agent Instructions)

You ran Phase C and got a report claiming:
- \(N(Z)=M\cdot D\) vanishes as \(\langle01\rangle^6\) (so \(\hat N=N/\langle01\rangle^6\) has deg ~24),
- cyclic brackets have double poles, non-cyclic simple poles,
- and **no** \(s_{ijk}\) channel poles (e.g. \(s_{123}\to0\) is finite/vanishing),
and concluded Phase C is “complete.”

**Do NOT accept “no channel poles” yet.** That conclusion is very likely a deformation/parameterization artifact.

This document clarifies what is solid, what must be re-tested, and the new north star.

---

## What seems solid (keep)

### A) \(\langle01\rangle\) vanishing order
If your clearing denominator satisfies:
- \(D(Z)\) is \(O(1)\) at \(\langle01\rangle\to0\),
- and the CHY amplitude suppresses strongly there,

then observing \(N(Z)=M\cdot D\sim \langle01\rangle^6\) is plausible.

✅ Keep the reduced object:
\[
\hat N \equiv \frac{N}{\langle01\rangle^6}.
\]

### B) “Cyclic double / non-cyclic simple” in *twistor coordinates*
It is plausible that in your chosen momentum-twistor parameterization, limits \(\langle i,i+1\rangle\to0\) look more singular because the map from twistors to spinors can include adjacent minors in denominators.

✅ Keep the divisor table as a *chart diagnostic*.

---

## What is likely wrong (must re-test)

### C) “No \(s_{ijk}\) poles” is **not trustworthy**
A tree-level gravity amplitude should factorize on physical multi-particle channels where:
\[
s_{123}=(p_1+p_2+p_3)^2 \to 0.
\]

If your probe says the amplitude is finite/vanishing there, the most likely explanations are:

1. Your deformation did **not actually approach a physical channel** (common in momentum twistors, since \(s_{ijk}\) is typically a ratio of bracket expressions).
2. Your deformation accidentally forces a different limit (collinear / chart singularity).
3. The object computed is not the full physical amplitude due to normalization/gauge/sector restrictions.

❌ Do **not** close Phase C until *channel factorization* is verified parameterization-invariantly.

---

## New North Star (Phase C completion criterion)

> **Phase C is not complete until you verify physical channel factorization.**

The agent’s job now:

1. Construct a deformation that guarantees **physical** \(s_{123}=\varepsilon\) (or solves for a parameter making it so),
2. Check that \(M_6\sim 1/\varepsilon\),
3. Match the residue to the product of two 4-pt amplitudes (sum over intermediate helicities):
\[
\mathrm{Res}_{s_{123}=0}M_6
\;=\;
\sum_{h=\pm} M_4(1,2,3,P^h)\,\frac{1}{s_{123}}\,M_4((-P)^{-h},4,5,6).
\]

If this passes robustly across random seeds and across all 3-particle channels, *then* Phase C is complete.

---

## Critical Fix: how to deform correctly

### Do NOT drive \(s_{123}\) by “making something small in \(Z\)” unless you verify it
Instead, compute \(s_{123}\) **directly from extracted spinors**:
- Build momenta \(p_i^{\alpha\dot\alpha}=\lambda_i^\alpha\tilde\lambda_i^{\dot\alpha}\).
- Compute \(s_{123}=(p_1+p_2+p_3)^2\) numerically/exactly at each \(\varepsilon\).

### Practical robust approach (recommended)
1. Pick a simple 1-parameter family \(Z(\varepsilon)\) (e.g. shift one twistor):
   - Example: \(Z_4 \mapsto Z_4 + \varepsilon Z_1\) (or similar).
2. For each \(\varepsilon\), compute:
   - \(s_{123}(\varepsilon)\),
   - \(M_6(\varepsilon)\).
3. Analyze the scaling of \(M_6\) vs \(s_{123}\):
   - Check whether \(s_{123}\cdot M_6\) approaches a **finite nonzero** limit as \(\varepsilon\to0\).

✅ **Pass condition:** \(s_{123}\cdot M_6\to c\ne 0\) and equals the factorization RHS.

---

## Concrete next task list (no wandering)

### Task C2-redux: Channel factorization oracle (start with \(s_{123}\))
For each 3-particle channel (start with 123):
1. Define a deformation family \(\varepsilon\mapsto Z(\varepsilon)\) that keeps adjacent \(\langle i,i+1\rangle\) generic (nonzero).
2. Compute \(s_{123}(\varepsilon)\) from spinors.
3. Compute \(M_6(\varepsilon)\).
4. Check:
   - scaling exponent of \(M_6\) vs \(s_{123}\),
   - \(\lim_{\varepsilon\to0}s_{123}(\varepsilon)M_6(\varepsilon)\) exists and is nonzero,
   - residue matches product of 4-pt oracles (Hodges 4pt).

### Required debugging outputs on failure
For any failed channel test, print:
- the measured scaling exponent \(\alpha\) in \(M_6\sim s_{123}^{\alpha}\),
- whether any adjacent brackets approached 0 (collinear/chart singularity),
- the deformation parameters (seeds + \(\varepsilon\) list),
- and compare LHS vs RHS normalization.

---

## Why the divisor table still matters
Keep the \(\langle ij\rangle\)-divisor table, but treat it as a **chart diagnostic**, not the physical geometry.

The physical positive-geometry target is tied to the **channel arrangement** \(\{s_S=0\}\), not raw \(\langle ij\rangle=0\) hyperplanes.

**Next mission:** translate from “bracket divisors” \(\to\) “channel divisors” and verify residues there.

---

## Stop Rules
1. Do **not** claim Phase C complete until channel factorization passes.
2. Do **not** trust any “\(s_{ijk}\to0\)” probe unless \(s_{ijk}(\varepsilon)\) is explicitly measured from spinors.
3. If factorization appears absent, assume the deformation missed the physical pole until proven otherwise.

