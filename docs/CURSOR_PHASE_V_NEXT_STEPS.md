# CURSOR — Phase V Next Steps (after Phase U): From "Map Sweep" to "Exact Jacobian & Residues"

> **Goal:** We have the correct geometry (Forest Polytope), the correct pushforward mechanism (Saddle Point), and the correct physics map (Spinor Edge Variables). The final piece is finding the **exact Jacobian factor** that normalizes the pushforward to equal the amplitude, and proving the **residue matching**.

---

## 0) Executive Summary of Phase U
*   **Success:** Implemented the "Saddle Point Pushforward" which correctly computes the canonical form of the polytope from the forest polynomial data (verified for n=4, 5).
*   **Success:** Classified all 26 facets of the n=6 Forest Polytope, matching them to physical subset constraints.
*   **Success:** Verified the "Core Identity" (Amplitude $\propto$ Forest Polynomial) is robust and reproducible.
*   **Gap:** The raw pushforward value $\Omega(X_{target})$ for n=6 is not yet equal to the amplitude magnitude (Ratio is $10^{-20}$ or erratic). This indicates we are missing a precise normalization factor (Jacobian of the map $z \to t$ or similar).

---

## 1) The Missing Factor: "Twisted" or "Weighted" Pushforward
The formula we implemented is:
\[ \Omega(X) = \sum \frac{1}{\det J} \]
But the stringy integral is actually:
\[ I = \int \frac{dx}{x} x^{\alpha X} p(x)^{-\alpha} \]
As $\alpha \to 0$, this becomes the canonical form. However, the *physical* object might be:
\[ M = \text{Prefactors} \times \Omega(X) \]
Or, the map itself involves coefficients that we set to 1.
The "Candidate A" map used $w_v \propto z^v$.
But in gravity, there are factors of $C_k^2$.
**Hypothesis:** The "weights" in the moment map should include the $C_k$ factors, or the Jacobian of the change of variables $z_{ij} \to \text{abstract parameters}$ contributes the missing factor.

## 2) Phase V Deliverables

### V1. The Normalization Search
*   **Action:** Modify `map_sweep.py` to test "Weighted Moment Maps".
    *   Instead of $c_v = 1$, set $c_v$ based on the graph structure (e.g., product of vertex weights?).
    *   Check if the ratio $M / \Omega$ becomes constant.
*   **Action:** Compute the "Hessian" of the log-polynomial explicitly and see if it matches the $C_k$ factors.

### V2. Residue Matching (The "Physics Proof")
*   **Action:** Pick a facet $F$ of the n=6 polytope (e.g., $x_{ij}=0$).
*   **Action:** Compute the residue of the saddle-pushforward on this facet.
    *   Mathematically: Res$(\Omega) = \Omega(\text{Facet})$.
*   **Action:** Verify this residue matches the $n=5$ amplitude (factorization).
*   **Code:** `src/posgeom/residue_checker.py`.

### V3. Full n=6 Numerical Match
*   **Goal:** Find the exact factor $K$ such that $M_{MHV} = K \cdot \Omega_{\text{Saddle}}$.
*   **Success Criteria:** CV of Ratio < 1% over 50 trials.

## 3) Immediate Next Commands
1.  `python src/scripts/map_sweep_weighted.py` (To be created)
2.  `python src/tests/test_residue_n5.py`



