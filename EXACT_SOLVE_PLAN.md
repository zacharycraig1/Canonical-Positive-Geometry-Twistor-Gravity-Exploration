# Reduced Numerator Solution

## Findings
The reduced numerator $P = N_H / \langle 0 1 \rangle^8$ was fit to a basis of 3 monomials:
1. $\langle 2 3 \rangle \langle 2 5 \rangle \langle 4 5 \rangle$ (Coeff ~ +85.15)
2. $\langle 2 4 \rangle \langle 2 5 \rangle \langle 3 5 \rangle$ (Coeff ~ -2.04)
3. $\langle 2 5 \rangle^2 \langle 3 4 \rangle$ (Coeff ~ -87.19)

## Issue
- The coefficients are **not simple integers** (e.g. 85.1486).
- The exact solve FAILED ("matrix equation has no solutions").
- This implies the 3-monomial basis is **insufficient** or the ansatz $N_H \propto \langle 0 1 \rangle^8$ is slightly wrong.
- Maybe there are other terms involving $0,1$ but with lower degree in $\langle 0 1 \rangle$?
- E.g. $\langle 0 1 \rangle^7 \langle 0 2 \rangle \dots$.
- But Weights were strictly [8,8,2,1,1,2].
- If term is $\langle 0 1 \rangle^7 \langle 0 2 \rangle \langle 1 5 \rangle$, degrees are:
  - 0: 7+1 = 8.
  - 1: 7+1 = 8.
  - 2: 1. (Need 2).
  - 5: 1. (Need 2).
  - Need more brackets for 2, 5.
- So we can have "leakage" from $\langle 0 1 \rangle$ to $\langle 0 X \rangle \langle 1 Y \rangle$.
- This suggests the basis is larger than 3 terms.
- **Hypothesis:** The "Reduced Numerator" is a polynomial in $\langle i j \rangle$ with weights $[8,8,2,1,1,2]$.
- We tested only terms with $\langle 0 1 \rangle^8$.
- We should allow $\langle 0 1 \rangle^7, \langle 0 1 \rangle^6$ etc.

## Final Plan
1.  Use the `reconstruct_numerator_v2.sage` solver (which generates ALL graphs for given degrees).
2.  It found 56 terms.
3.  The fit was exact (residuals empty).
4.  The coefficients were "messy".
5.  **Action:** Use Sage's `matrix.solve_right` on the 56-term system to get **exact rational coefficients**.
6.  The result will be the **exact analytic formula** for the Hodges Numerator.
7.  Once we have the rational coefficients, we can try to simplify them using Schouten identities (or just report the rational form).








