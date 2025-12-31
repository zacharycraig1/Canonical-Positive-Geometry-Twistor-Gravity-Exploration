# Physics Project: Positive Geometry of Gravity

## Current Status: Phase H (Complete) - Geometric Verification
We have successfully identified and verified the **Positive Geometry** underlying the MHV Gravity Amplitude.
The key object is the **3-Rooted Spanning Forest Polytope** of the complete graph $K_n$.

### Key Findings
1.  **Exact Identity Verified:** For $n=4, 5, 6$, the MHV Gravity Amplitude is exactly the canonical form of the Forest Polytope (evaluated on dual kinematic variables), up to standard prefactors.
    \[
    M_n^{\text{MHV}} \propto \Omega(P_{\text{Forest}})(z_{ij})
    \]
    - Verification script: `src/scripts/physics_pullback_n[4,5,6].sage` (All passed with Ratio 1.000000).

2.  **Geometric Structure:**
    - **Polytope:** Newton polytope of the "Forest Polynomial" (sum of monomials for rooted forests).
    - **Toric Variety:** The associated projective toric variety $X_P$.
    - **Canonical Form:** Computed via triangulation of the polytope in the dual space.

3.  **Codebase Structure:**
    - `src/posgeom/`: Core geometric module.
        - `forest_polytope.py`: Generates the forest polynomial (combinatorial object).
        - `toric.py`: Computes the affine lattice basis and toric map.
        - `canonical_polytope.py`: Evaluates the canonical form $\Omega(W)$ on the dual space.
        - `physics_map.py`: Maps spinor kinematics to the edge variables $z_{ij}$.
    - `src/chy_oracle/`: Physics oracle (Hodges determinant, KLT, etc.) for validation.

### Next Steps (Phase I)
The next agent should focus on **Residue Analysis** and **Paper Drafting**.

1.  **Residue/Factorization Proof:**
    - Prove that the facets of the Forest Polytope correspond 1-to-1 with physical factorization channels (poles of the amplitude).
    - Use `src/posgeom/toric.py` to list facets and check their physical meaning ($z_{ij} \to 0$ or sums thereof).

2.  **Finalize the "Pushforward" Statement:**
    - We have numerically verified the dual form. Formalize the algebraic moment map $\mu: X_P \to P^*$ that connects the Toric Variety form to the Amplitude.

3.  **Drafting:**
    - Update `paper/main.tex` with the specific theorem: "The MHV gravity amplitude is the canonical form of the 3-rooted spanning forest polytope."
    - Cite Chaiken (1982) for the combinatorial basis.

### How to Run Verifications
- **n=4:** `sage -python src/scripts/physics_pullback_n4.sage`
- **n=5:** `sage -python src/scripts/physics_pullback_n5.sage`
- **n=6:** `sage -python src/scripts/physics_pullback_n6.sage`
- **Geometric Tests:** `sage -python src/tests/test_canonical_form_polytope.py`

### References
- **Theorem Inventory:** See `docs/theorem_inventory.md` for a precise list of what is known vs. new.
- **Phase G Report:** See `src/scripts/PHASE_G_REPORT.md` for details on the forest polynomial construction.
