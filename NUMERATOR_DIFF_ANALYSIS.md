# Numerator Difference Analysis

## Findings
- **Numerator of Hodges ($N_H$):** Confirmed to be a polynomial of degree 34 (Weight 22/34?).
- **Difference ($N_H - N_K$):** Failed to extrapolate as a polynomial.
  - Extrapolation check had relative error ~0.3% ($3.26e10$ vs $3.25e10$).
  - Degree 39 fit (on 40 points).
  - This suggests the difference is **NOT** a polynomial, or has very high degree.
  - **Conclusion:** $N_{KLT} = M_{KLT} \times D_{cyclic}$ is **NOT** a polynomial.
  - This means $M_{KLT}$ has poles that are NOT in $D_{cyclic}$.
  - $D_{cyclic}$ contains only $\langle i i+1 \rangle$.
  - $M_{KLT}$ likely contains other poles (e.g. $\langle i j \rangle$ or $s_{ij}$).

## Implication
The standard KLT formula introduces spurious poles (likely $s_{ij}$ poles from the kernel or partial amplitudes).
In the correct amplitude, these spurious poles must cancel.
Since $M_{Hodges}$ has only physical poles (confirmed by $N_H$ being polynomial over physical denominator), the KLT sum as implemented fails to cancel these poles.

## Next Steps
1. **Identify Spurious Poles in KLT:** Check which $s_{ij}$ poles remain in $M_{KLT}$.
2. **Fix KLT Basis:** The failure to cancel poles suggests the BCFW/KLT basis or the Kernel is wrong.
3. **Use Hodges N as Target:** Since we know $N_H$ is the correct polynomial, we should try to expand it directly in a basis of local functions.







