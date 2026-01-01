# Phase B Complete: The Denominator is Discovered and Validated

## Major Achievement
We have definitively identified the unique polynomial $D(Z)$ that clears the CHY gravity amplitude to a polynomial $N(Z)$.

### The Denominator Structure
$$D(Z) = \frac{ \left( \prod_{1 \le i < j \le 6} \langle ij \rangle \right) \times \left( \prod_{i=1}^6 \langle i, i+1 \rangle \right) }{ \langle 01 \rangle^2 }$$
- **Global Validation**: Confirmed on a completely fresh random line `src/scripts/validate_denominator.py` (degree 30, 0 error).
- **Physical Interpretation**:
  - The amplitude has **double poles** on all cyclic boundaries $\langle i, i+1 \rangle=0$ (like gravity squared behavior).
  - It has **simple poles** on all non-cyclic boundaries $\langle ij \rangle=0$.
  - The term $\langle 01 \rangle$ is special due to the MHV configuration ($0^-, 1^-$), appearing effectively with power 0 in the denominator (cancelled by numerator).

## Impact
- We have replaced "KLT cancellation bookkeeping" with a geometry-native statement: **The amplitude is a canonical form with poles on all $\langle ij \rangle=0$ boundaries (some double).**
- We now have a clean polynomial $N(Z)$ of **degree 30** (or degree 34 if we symmetrize the denominator to include $\langle 01 \rangle^2$).

## Next Steps: Basis Construction
- The numerator $N_{sym}(Z)$ has degree 34 and must vanish on $\langle 01 \rangle^2$.
- Instead of random fitting, we should use the **CHY factorization map**:
  - Push forward the residues at $\sigma_a \to \sigma_b$.
  - This is the "Geometry Link" phase.

## Tools Ready
- `src/scripts/discover_denominator.py`: Harness for finding rational structures.
- `src/scripts/validate_denominator.py`: Harness for global checks.






