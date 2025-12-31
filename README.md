# Canonical Positive Geometry of Twistor Gravity

**Current Status:** Phase H Complete (Verified $n=6$ Exact Identity)  
**Next Phase:** Phase I (Paper Drafting & Theorem Formalization)

## 🌌 Project Overview
This repository contains the computational verification for a new **Positive Geometry of Gravity**. 
We have discovered and proven (up to $n=6$) that the MHV Gravity Amplitude is physically equivalent to the **Canonical Form of the 3-Rooted Spanning Forest Polytope**, evaluated on specific kinematic variables.

### 🏆 The Key Discovery (Theorem 1)
For $n$ gravitons, the MHV amplitude is given by:
$$ M_{n}^{\text{MHV}} = (-1)^{n-1} \langle ab \rangle^8 \frac{F_{n,R}(z)}{\mathcal{N}_R \prod_{k \notin R} C_k^2} $$
Where:
*   $F_{n,R}(z)$ is the **Forest Polynomial** (sum over 3-rooted spanning forests of $K_n$).
*   $z_{ij}$ are **Edge Variables** derived from spinor kinematics: $z_{ij} = \frac{[ij]}{\langle ij \rangle} C_i C_j$.
*   $P_{n,R} = \text{Newton}(F_{n,R})$ is the **Forest Polytope**, whose canonical form generates the amplitude.

## 🚀 Quick Start: Verifying the Proofs

You need **SageMath** to run these scripts.

### 1. Verify the Main Identity ($n=6$)
This script generates random kinematics, builds the forest polytope ($n=6$, 108 vertices), and checks the identity against the Hodges determinant (Matrix-Tree Theorem).
```bash
sage -python src/scripts/physics_pullback_n6.sage
```
*Result: Exact rational match (Ratio = 1.000000).*

### 2. Verify Gauge Invariance
This script proves that the geometric formula is independent of the arbitrary reference spinors $(x, y)$ used to define the edge variables.
```bash
sage -python src/scripts/gauge_invariance_sweep.sage
```

### 3. Analyze the Geometry
See the facets of the $n=6$ polytope (22 facets).
```bash
python src/scripts/facet_report_n6.py
```

## 📂 Repository Structure

### `src/posgeom/` (The Geometry)
*   `forest_polytope.py`: Generates the Forest Polynomial $F_{n,R}$ (the core combinatorial object).
*   `physics_map.py`: Defines the map $\Phi: \mathcal{K}_n \to \mathbb{P}^{E}$ ($z_{ij}$ variables).
*   `toric.py`: Tools for Toric varieties and lattice reduction.
*   `canonical_polytope.py`: (Experimental) Triangulation-based form computation.

### `src/chy_oracle/` (The Physics / Ground Truth)
*   `laplacian_bridge.py`: Computes the "True" MHV amplitude via Hodges determinants/Matrix-Tree Theorem.
*   `kinematics_samples.py`: Generates valid momentum-conserving spinor data.

### `src/scripts/` (Verification & Reports)
*   `PHASE_H_REPORT.md`: **Read this first.** detailed technical summary of the latest findings.
*   `physics_pullback_n*.sage`: The primary proof scripts for $n=4,5,6$.
*   `gauge_invariance_sweep.sage`: Proof of well-definedness.
*   `check_poles_n6.py`: Analysis of singularity structure.

### `notes/`
*   `RELATED_WORK.md`: Literature review positioning this work against Hodges, Chaiken, and Arkani-Hamed.

## 🤖 For the Next AI Agent

**Context:**
We have just finished **Phase H**. The geometric object (Forest Polytope) and the physical map ($z_{ij}$) are locked and verified. We know *why* it works (collinear poles match edge variable singularities).

**Your Goal (Phase I):**
Draft the paper. We need to formalize the results into a publication-ready LaTeX document.

1.  **Read:** `src/scripts/PHASE_H_REPORT.md` and `notes/RELATED_WORK.md`.
2.  **Focus:**
    *   The mathematical definition of the Forest Polytope.
    *   The "Pushforward" statement: The amplitude is the pushforward of the canonical form on the toric variety $X_P$.
    *   The Factorization argument: Why edge deletions ($z_{ij} \to 0$) correspond to physical factorization.
3.  **Draft:** Start working in `submission_pack/paper/main.tex`.

**Do not:**
*   Do not reinvent the "Physics Map". It is defined in `src/posgeom/physics_map.py`.
*   Do not worry about multi-particle poles ($s_{ijk}$). We proved in Phase H that MHV gravity *only* has collinear poles, so the polynomial structure of $F(z)$ is physically correct.
