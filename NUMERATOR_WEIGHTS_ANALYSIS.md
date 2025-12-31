# Hodges Numerator Weights Analysis

## Findings
The weights of the Hodges Numerator $N_H(Z)$ with respect to each particle $Z_i$ are:
- $Z_0$: Weight 8
- $Z_1$: Weight 8
- $Z_2$: Weight 2
- $Z_3$: Weight 1
- $Z_4$: Weight 1
- $Z_5$: Weight 2

**Total Weight:** 22.

## Interpretation
- **Asymmetry:** The weights are highly asymmetric.
  - Particles 0,1 (Legs 1,2) have Weight 8.
  - Particles 2,5 (Legs 3,6) have Weight 2.
  - Particles 3,4 (Legs 4,5) have Weight 1.
- **Reason:** This asymmetry arises because the Hodges formula (as implemented) depends on a **deletion set** (here $\{0,1,2\}$? or $\{0,4,5\}$?) and explicitly breaks symmetry (restored only when dividing by normalization factors).
- Wait, the Hodges amplitude *should* be symmetric.
- But $N_H = M_H \times D_{cyclic}$. $D_{cyclic}$ is symmetric.
- So $N_H$ *should* be symmetric.
- Why are the weights asymmetric?
  - **Crucial:** The `hodges_6pt_mhv` function uses a *fixed deletion set* (default `[0,1,2]`) unless specified.
  - In `check_weights.sage`, I called `hodges_6pt_mhv(tw)`. It uses default deletion.
  - The resulting $M_H$ is permutation invariant (proven earlier).
  - So $N_H$ should be symmetric.
  - **Unless** the implementation of `hodges_6pt_mhv` introduces asymmetry via reference spinors or normalization?
  - But `check_proper_symmetry` showed Ratio=1.0000. So $M_H$ IS symmetric.
  - If $M_H$ is symmetric and $D_{cyclic}$ is symmetric, $N_H$ MUST be symmetric.
  - So weights MUST be equal. $22/6 \approx 3.66$. Not integer.
  - **Contradiction:** Weights are integers 8, 8, 2, 1, 1, 2.
  - This implies $M_H$ is **NOT** symmetric under independent scaling of $Z_i$?
  - $Z_i \to c Z_i$ is a Little Group scaling.
  - For Graviton MHV, Little Group weights are uniform (spin 2).
  - $M(c \lambda, c^{-1} \tilde\lambda) = c^{-2h} M$.
  - In momentum twistors, $Z \to c Z$ is not just Little Group. It scales the geometry.
  - **Wait:** In twistors, $M$ has weight 0 in $Z_i$?
  - No, dual conformal weight.
  - Graviton amplitude has weights.
  - The weights (8,8,2,1,1,2) sum to 22.
  - This asymmetry suggests the "Hodges Amplitude" I am computing might be in a specific "gauge" or "frame" defined by the reference spinors or deletion set?
  - But $M$ is unique.
  - **Explanation:** The factor $\langle 0 1 \rangle^8$ in `hodges_6pt_mhv`!
  - `helicity_factor = twistor.get_angle(0, 1) ** 8`.
  - This factor gives Weight 8 to $Z_0, Z_1$.
  - This breaks symmetry explicitly!
  - Why is it there?
  - "For MHV with 0,1 negative helicity".
  - Ah! MHV amplitudes are only symmetric up to helicity factors.
  - $M_{MHV} = \frac{\langle i j \rangle^8}{\dots}$? No, that's N=8 Supergravity?
  - Pure Gravity MHV: $M = \langle i j \rangle^8 \dots$
  - The indices $i,j$ are the negative helicity gravitons.
  - So the numerator depends on which particles are negative helicity.
  - This explains the asymmetry!
  - Particles 0 and 1 are the negative helicity gravitons.
  - So they carry extra weight.

## Implication for Polynomial Basis
To reconstruct $N_H$, we need a basis that respects these specific weights:
- Deg(0)=8, Deg(1)=8, Deg(2)=2, Deg(3)=1, Deg(4)=1, Deg(5)=2.
- Variables: 4-brackets $\langle abcd \rangle$.
- Each bracket contributes 1 to four indices.
- We need to form a monomial $\prod \langle a b c d \rangle^{n_{abcd}}$ such that the degree in $Z_i$ matches the target.
- This is a Diophantine problem.
- Also need to consider $\langle i j \rangle$ (Weight 2 in Z) if we work with angle brackets.
- The numerator is likely expressible in angle brackets (Spinor Helicity).
- Or 4-brackets (Momentum Twistors).
- Since $D_{cyclic}$ uses angle brackets, $N_H$ likely uses angle brackets.
- Angle bracket basis: Monomials in $\langle i j \rangle$.
- Constraints:
  - $\sum_{j \ne i} n_{ij} = w_i$.
  - $w = [8, 8, 2, 1, 1, 2]$.
  - $n_{ij} \ge 0$.
- Let's construct this basis.

## Plan
1.  **Generate Basis:** Find all graphs (multigraphs) with degree sequence $[8, 8, 2, 1, 1, 2]$.
    - Nodes 0..5. Edges $(i,j)$ correspond to $\langle i j \rangle$.
    - Number of edges = Total Weight / 2 = 22 / 2 = 11.
    - Find all multigraphs with 11 edges and these degrees.
2.  **Fit:** Use this basis to fit $N_H$.
3.  **Result:** The coefficients will give the explicit formula.




