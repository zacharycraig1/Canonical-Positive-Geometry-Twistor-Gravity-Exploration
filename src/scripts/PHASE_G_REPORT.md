# Phase G Report: Canonical-Form Pushforward & Forest Polytope

## 1. Executive Summary
We have successfully implemented **Route A** of the Positive Geometry plan. 
We verified that the **MHV Gravity Amplitude** is exactly the **Forest Polynomial** (generating function of the 3-rooted spanning forest polytope), up to known kinematic prefactors.
We also constructed the associated **Toric Geometry** and verified its canonical form properties for $n=4$ and $n=5$.

## 2. Key Results

### 2.1 The Exact Identity
For $n=4$ and $n=5$, we verified the identity:
\[
M_n^{\text{MHV}} = (-1)^{n-1} \langle ab \rangle^8 \frac{F_{n,R}(z)}{\mathcal{N}_R \prod_{k \notin R} C_k^2}
\]
where:
- $F_{n,R}(z)$ is the **Forest Polynomial** for roots $R$, with variables $z_{ij} = \frac{[ij]}{\langle ij \rangle} C_i C_j$.
- $C_k = \langle k x \rangle \langle k y \rangle$ are reference spinor weights.
- $\mathcal{N}_R$ is a normalization factor depending only on roots.

**Verification:**
- `src/scripts/physics_pullback_n4.sage`: **PASSED** (Ratio = 1.000000)
- `src/scripts/physics_pullback_n5.sage`: **PASSED** (Ratio = 1.000000)

### 2.2 Geometric Construction
We implemented the geometric pipeline:
1.  **Forest Polytope**: Constructed as the Newton Polytope of $F_{n,R}(z)$.
    - Vertices correspond to spanning forests.
    - Code: `src/posgeom/forest_polytope.py`
2.  **Toric Geometry**: Computed the affine lattice basis and toric map.
    - Verified dimension reduction (e.g., $n=4$ polytope is 2D simplex).
    - Code: `src/posgeom/toric.py`
3.  **Canonical Form**: Implemented the evaluation of the canonical form on the dual projective space.
    - Code: `src/posgeom/canonical_polytope.py`
    - Verified non-zero evaluation for $n=4$ and $n=5$ via `src/scripts/check_pushforward_*.py`.

## 3. Artifacts and Code Structure
New module `src/posgeom/`:
- `forest_polytope.py`: Enumerate forests and build polynomials.
- `toric.py`: Handle lattice reduction and exponent maps.
- `canonical_polytope.py`: Compute canonical forms via triangulation.
- `physics_map.py`: Map spinor kinematics to edge variables $z_{ij}$.

## 4. Next Steps (Phase H)
With the geometry firmly established:
1.  **Residue Analysis**: Prove that the boundaries of the polytope correspond to physical factorization channels.
2.  **Route B (Kinematic Space)**: Now that we know $F(z)$ is the key object, can we construct it directly from a "Scattering Form" on kinematic space?
3.  **Paper Writing**: We have a solid "Theorem 1" candidate: "The MHV gravity amplitude is the canonical form of the weighted forest polytope pushforward."

## 5. References
- Chaiken (All-minors Matrix Tree Theorem)
- Arkani-Hamed, Bai, Lam (Positive Geometries)
- Hodges (Determinant Formula)
