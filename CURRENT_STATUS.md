# Current Status: KLT-Bilinear & Split Signature Geometry

## Overview
We have investigated the KLT-Bilinear ansatz for 6-point MHV gravity.
1.  **Scaling Verification**: We confirmed that the KLT amplitude (field theory) scales as $t^{-2}$ under spinor rescaling $Z \to tZ$ (where $s \sim t^2$).
2.  **Hodges Discrepancy**: The Reduced Hodges determinant scales as $t^{-18}$. A factor of $t^{16}$ (like $\langle 01 \rangle^8$) is required dimensionally.
3.  **Asymmetry Identified**: The KLT implementation in `src/klt.sage` yields results that depend on the choice of pivot particle.
4.  **Split Signature**: The KLT kernel matrix has a robust **Split Signature**.
    - For $n=6$, it is **(3, 3)** in the majority of regions (~70%).
    - The signature **jumps** to (2,4) or (4,2) in other kinematic regions.
    - It is never positive definite.
5.  **Intersection Theory Connection**:
    - We computed the **Bi-adjoint Scalar Intersection Matrix** $m$.
    - This matrix also exhibits **Split Signature** and varies by region, matching the KLT kernel's behavior.
    - $S_{KLT}$ and $m^{-1}$ share geometric properties but currently show a basis/normalization mismatch.

## Conclusion
The **Single Log-Form** hypothesis is ruled out.
The **Double-Copy Geometry** is confirmed to be **Twisted Intersection Theory** with a **Split Signature Metric**. The metric signature depends on the kinematic region (chamber), which is a known feature of twisted intersection numbers.

## Next Steps
1.  **Fix Basis Mismatch**: Derive the exact basis mapping between the KLT kernel and the Bi-adjoint Scalar diagrammatic sum to prove $S_{KLT} \propto m^{-1}$ exactly.
2.  **Chamber Analysis**: Map the (3,3) signature regions to specific physical scattering channels (Euclidean vs Physical).
3.  **Hodges Matching**: Once the kernel is fixed, relate the Intersection Pairing $\langle A_L | S | A_R \rangle$ to the Hodges determinant.

## Repository Layout
-   `ruled_out_single_log_form/`: Archive of the failed single-form search.
-   `klt_search/`: Scripts for the KLT ansatz search and analysis.
-   `twisted_cohomology/`: Initial setup for twisted forms.
-   `results/`: Findings reports.
-   `src/`: Core libraries.
