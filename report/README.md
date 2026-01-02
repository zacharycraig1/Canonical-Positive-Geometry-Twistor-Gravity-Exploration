# Physics Project Report: Positive Geometry of MHV Gravity

This folder contains the complete verification suite and results for the project identifying the **Weighted Laplacian / 3-Rooted Forest** geometry underlying the MHV Gravity Amplitude.

## 📂 Directory Structure

### `src/notes/` (Key Findings)
- `RESULTS_MASTER.md`: **Start Here.** Summary of the core Master Identity and geometric interpretation.
- `PHASE_F_THEOREM.md`: Detailed mathematical statement of the theorem.
- `PHASE_E_REPORT.md`: Report on the exact algebraic bridge and pole analysis.

### `src/chy_oracle/` (Core Libraries)
- `amplitude_spinor.py`: Generalized $n$-point Hodges MHV amplitude implementation.
- `laplacian_bridge.py`: The bridge code reconstructing the amplitude from the Weighted Laplacian.
- `forest_sum.py`: Combinatorial implementation of the All-Minors Matrix-Tree Theorem (Forest Sum).
- `matrix_tree.py`: Matrix-Tree Theorem and Weighted Laplacian construction.
- `kinematics_samples.py`: Momentum-conserving kinematics generator (using Twistors).

### `src/scripts/` (Verification Scripts)
- **Generalization ($n=7$):**
    - `phaseF3_n7_verification.py`: Verifies the Master Identity for $n=7$ (Passes).
    - `phaseF4_n7_valuations.py`: Analyzes pole structure for $n=7$.
- **Robustness:**
    - `phaseF1_reference_independence.py`: Verifies gauge invariance (independence of reference spinors).
- **Geometry:**
    - `phaseF2b_forest_newton_polytope.py`: Computes the Newton Polytope of the Forest Sum.
- **Foundational ($n=6$):**
    - `phaseE3_reconstruct_M_from_tree.py`: Verifies the Master Identity for $n=6$.
    - `phaseE3_chy_to_tree_bridge.py`: Verifies the link to CHY Jacobian.
    - `phaseE1_tree_vs_detprime.py`: Pins the exact normalization factors.
    - `phaseE2_spinor_valuations.py`: Detailed pole analysis for $n=6$.

## 🚀 How to Run

You need a Python environment with **SageMath** installed.

To verify the main result ($n=7$ extension):
```bash
sage -python src/scripts/phaseF3_n7_verification.py
```

To verify gauge invariance:
```bash
sage -python src/scripts/phaseF1_reference_independence.py
```

To verify the Forest Sum theorem ($n=6$):
```bash
sage -python src/chy_oracle/forest_sum.py
```







