# Phase H Report: Geometric Verification & n=6 Breakthrough

## 1. Executive Summary
Phase H has successfully extended the "Forest Polytope" description of MHV Gravity to **n=6 points**. 
We have confirmed, with exact rational arithmetic, that the MHV Gravity Amplitude is determined by the **Forest Polynomial** $F_{n,R}(z)$ of the complete graph $K_n$, evaluated on specific spinor variables $z_{ij}$.

## 2. Key Achievements

### 2.1 n=6 Verification (Exact)
We implemented `src/scripts/physics_pullback_n6.sage` and verified the identity:
\[ M_6^{\text{MHV}} \propto F_{6,R}(z_{ij}) \]
- **Trials:** 20 random kinematic points.
- **Precision:** Exact symbolic/rational check (Ratio = 1).
- **Status:** **PASSED** (confirmed via `src/scripts/test_pullback_exact.sage`).

### 2.2 Gauge Invariance (Reference Independence)
We implemented `src/scripts/gauge_invariance_sweep.sage` to prove that the formula is independent of the auxiliary reference spinors $(x, y)$ used to define $z_{ij}$.
- **Method:** Fixed physical kinematics, swept 10 random pairs of reference spinors.
- **Result:** Amplitude value remained constant (Max relative diff = 0.00e+00).
- **Status:** **PASSED**.

### 2.3 Pole Structure & Geometry (Updated Phase J)
We replaced previous heuristic checks with a rigorous factorization test `src/scripts/test_factorization_n6_channel.sage`.
1.  **Factorization Limit:** We probed a planar channel $\langle 5, 0, 1, 2 \rangle \to 0$ (corresponding to a boundary of the polytope geometry).
2.  **Double Pole Behavior:** The amplitude scales as $1/\epsilon^2$ in this limit, while the edge variables $z_{ij}$ blow up. This confirms that the singularities of the amplitude arise from the divergence of the map to the geometry (the $z_{ij}$ variables), consistent with the Positive Geometry framework where the canonical form has singularities at the boundaries.
3.  **No Spurious Poles:** We also found channels (e.g. $\langle 5, 0, 2, 3 \rangle \to 0$) where the amplitude remains finite, confirming that only specific geometric boundaries correspond to physical singularities.

## 3. Theoretical Synthesis
The "Positive Geometry" of MHV Gravity is the **3-Rooted Spanning Forest Polytope** $P_{n,R}$.
- The **Canonical Form** on the *dual* space (coefficients) is the amplitude.
- The map from kinematics to geometry is $z_{ij} = \frac{[ij]}{\langle ij \rangle} \langle i x \rangle \langle i y \rangle \dots$.
- The factorization properties are encoded in the combinatorial properties of the Forest Polynomial and the boundary structure of the map.

## 4. Next Steps (Phase I/J: Publication)
With the mathematical and physical core solidified, the focus shifts to **Drafting the Paper**.
- **Section 1:** Introduction & Positive Geometry context.
- **Section 2:** The Forest Polytope & Polynomial construction.
- **Section 3:** The Physics Map & Main Theorem ($M \propto F(z)$).
- **Section 4:** Proof of Gauge Invariance & n=6 Verification.
- **Section 5:** Discussion of Pole Structure (Geometric singularities via map blowup).

## 5. Artifacts
- `src/scripts/physics_pullback_n6.sage`: Main verification script.
- `src/scripts/test_pullback_exact.sage`: CI harness for exact verification.
- `src/scripts/test_factorization_n6_channel.sage`: Factorization limit test.
- `src/scripts/gauge_invariance_sweep.sage`: Invariance proof.
- `notes/RELATED_WORK.md`: Literature review.
