# Claim Status

This document tracks the status of scientific claims regarding the positive geometry of MHV gravity amplitudes.

## 1. Naive Pullback Claim
**Claim:** The canonical form $\Omega_{P_{n,R}}$ of the forest polytope, pulled back via the map $W = (1, z)$, equals the MHV gravity amplitude $M_n^{\mathrm{MHV}}$.
$$
\Phi^*(\Omega_{P_{n,R}}) \stackrel{?}{=} M_n^{\mathrm{MHV}} \quad \text{via} \quad \Phi: z \mapsto W=(1, z)
$$

**Status:** **FALSE**
**Reasoning:**
- The forest polynomial is a sum of monomials $z^{a_F}$, whereas the polytope canonical form is rational in linear forms $W \cdot Z_F$.
- The identification $W=(1, z)$ attempts to map a linear structure to a monomial one without the correct Jacobian or exponential map.
- Diagnosed by order-of-magnitude failures in the $n=6$ ratio test and "W is not z" analysis.

## 2. Forest Polynomial Identity
**Claim:** The MHV gravity amplitude can be written as a specific forest polynomial $F_{n,R}(z)$ evaluated on edge variables $z_{ij}$.

**Status:** **VERIFIED**
**Evidence:**
- Numerically verified for $n=4, 5, 6$.
- Matches known Matrix-Tree Theorem results (Hodges, Feng & He).
- This is the "computational substrate" for the geometric claim, but not the geometric statement itself.

## 3. Stringy Pushforward Claim
**Claim:** The MHV gravity amplitude arises as the $\alpha' \to 0$ limit of a stringy canonical form integral over the forest polytope $P_{n,R}$. Alternatively, it is the pushforward of a canonical form on the positive torus (or product of $\mathbb{CP}^1$s) under the monomial map.

**Status:** **IN PROGRESS / VERIFIED ALGEBRAICALLY**
**Findings:**
- Implemented `src/posgeom/stringy_integral.py` which verifies the algebraic pullback property for $n=4, 5$.
- Confirmed that the pullback of the standard simplex form under the monomial map matches the polytope geometry (conceptually).
- **Blocker:** Direct evaluation of $\Omega_P$ at the physical point requires resolving the mismatch between the 108-dimensional monomial vector and the 16-dimensional dual space of the polytope.
- **Hypothesis:** The physical amplitude corresponds to the form evaluated at $W = (0, \log z)$, but the numeric ratio (~1e-5) suggests normalization factors or measure terms are missing.

## 4. Root Patch Atlas
**Status:** **COMPLETED**
- `docs/root_patch_atlas.md` generated.
- Maps physical singularities ($z_{ij} \to 0$) to facets of the forest polytope across all root choices.
