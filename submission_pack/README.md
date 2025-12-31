# n-Point MHV Gravity: Positive Weighted Laplacian Geometry

**We have identified and proven the geometric object underlying the $n$-point MHV Gravity Amplitude: the 3-Rooted Spanning Forest Polytope of $K_n$.**

The amplitude is computed exactly by the principal minor of a **Weighted Laplacian Matrix** (corank 3), which expands combinatorially into a sum over **3-rooted spanning forests**.

This construction has been verified for $n=6$ and $n=7$ and proves that gravity amplitudes are "Positive Geometries".

---

## 🚀 Reproduction (Publication Pack)

A complete verification suite is available in the **`repro/`** directory.

To run all checks (Verification, Valuations, Combinatorics, Polytopes):

```bash
# Linux / Mac
./repro/run_all.sh

# Windows (PowerShell)
.\repro\run_all.ps1
```

### Verification Manifest
1.  **Reference Independence:** `src/scripts/phaseF1_reference_independence.py`
2.  **Deletion Independence:** `src/scripts/phaseF2_deletion_set_independence.py`
3.  **Pole Order Audit:** `src/scripts/phaseF3_exact_pole_orders.py`
4.  **Forest Combinatorics:** `src/scripts/phaseF4_all_minors_forest_expansion.py`
5.  **Polytope Geometry:** `src/scripts/phaseF5_newton_polytopes.py`
6.  **n=7 Generalization:** `src/scripts/phaseF6_n7_verification.py`

---

## 📖 Key Formula

$$ \mathcal{M}_n = (-1)^{n-1} \langle ab \rangle^8 \cdot \frac{\det(\tilde{L}^{(R)})}{\prod_{k \notin R} C_k^2 |\mathcal{N}(R)|^2} $$

Where $\tilde{L}$ is the weighted Laplacian with weights $w_{ij} = [ij]/\langle ij \rangle$ and vertex factors $C_i = \langle i x \rangle \langle i y \rangle$.

---

## 📂 Structure

- `src/`: Core python/sage library (`chy_oracle`).
- `results/`: JSON logs of verification runs.
- `paper/`: LaTeX skeleton of the publication.
- `repro/`: One-click reproduction scripts.

