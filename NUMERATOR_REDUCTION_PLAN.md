# Numerator Reconstruction Analysis

## Findings
The least-squares fit produced a formula with 56 terms (all basis elements used), with **non-integer** coefficients (e.g., -1660.1466, +3440.6830).
This suggests that the basis I constructed (via random search or Diophantine solver) is **not the correct minimal basis** or that the coefficients are rational/large integers.

However, one term stands out:
$$ +8062.2293 \times \langle 0 1 \rangle^8 \langle 2 5 \rangle^2 \langle 3 4 \rangle $$
This corresponds to my "Single Monomial Hypothesis" (Weights 8,8,2,1,1,2).
But there are many other large terms, e.g.
$$ +6739.0549 \times \langle 0 1 \rangle^7 \langle 0 5 \rangle \langle 1 2 \rangle \langle 2 3 \rangle \langle 4 5 \rangle $$
(Need to parse the string better).

## Diagnosis
- The "Weights" analysis was correct (8,8,2,1,1,2).
- The "Basis" construction was exhaustive for the given weights?
- `reconstruct_numerator_v2.sage` found 56 monomials.
- The residuals were empty `Residual: []`. This means exact fit on the sample points (since I used num_basis + 10 points).
- **BUT** the coefficients are not simple integers.
- This implies either:
  1. The "Angle Bracket Basis" is not the natural basis. Maybe 4-brackets?
  2. The coefficients are fractions with large denominators?
  3. There is linear dependence in the basis (Angle bracket identities, Schouten identities).
- Since we are working with angle brackets $\langle i j \rangle$ for $n=6$, there are Schouten identities (3-term relations) and cyclic identities.
- The space of weight-22 monomials is overcomplete.
- The solver found *a* solution, but not a canonical one.

## Refinement
The fact that coefficients are large and messy strongly suggests we should use **Schouten Identities** to reduce the basis to a minimal one.
Alternatively, we can try to find a "nice" representation.

Notice the structure of the "dominant" terms:
- Many terms have $\langle 0 1 \rangle^7$ or $\langle 0 1 \rangle^6$.
- The hypothesis term $\langle 0 1 \rangle^8 \langle 2 5 \rangle^2 \langle 3 4 \rangle$ has coefficient ~8062.
- The amplitude values were around $\sim 1000$.
- So coefficients are order 1.
- Why 8062?
- Maybe normalization is off?
- `hodges` result has `<0 1>^8` factor.
- Maybe we should divide N by `<0 1>^8`?
- Then we look for a polynomial of weight [0,0,2,1,1,2] + [0,0,0,0,0,0] = [0,0,2,1,1,2].
- Total Weight 6.
- This is MUCH simpler!
- Let's try to fit $N' = N / \langle 0 1 \rangle^8$.
- Target Weights: 0, 0, 2, 1, 1, 2.
- Monomials must NOT contain 0 or 1?
- No, weights are degrees in $Z_i$.
- $\langle 0 1 \rangle^8$ contributes 8 to $Z_0$ and 8 to $Z_1$.
- If we divide by it, the remaining factor $P$ must have:
  - Deg(0) = 8 - 8 = 0.
  - Deg(1) = 8 - 8 = 0.
  - Deg(2) = 2.
  - Deg(3) = 1.
  - Deg(4) = 1.
  - Deg(5) = 2.
- This means $P$ depends ONLY on $Z_2, Z_3, Z_4, Z_5$!
- **Huge Simplification!**
- The polynomial $P$ is a function of only 4 twistors.
- And it has very low degree (Total 6).
- Basis: Monomials in $\langle 2 3 \rangle, \langle 2 4 \rangle, \langle 2 5 \rangle, \langle 3 4 \rangle, \langle 3 5 \rangle, \langle 4 5 \rangle$.
- Constraints:
  - Deg(2)=2
  - Deg(3)=1
  - Deg(4)=1
  - Deg(5)=2
- Possible edges (on nodes 2,3,4,5):
  - (3,4) is likely (uses deg 1,1).
  - Remaining: 2 needs 2, 5 needs 2.
  - (2,5)^2 uses them.
  - So $\langle 3 4 \rangle \langle 2 5 \rangle^2$ is a candidate.
  - Other options:
    - (2,3), (4,5) ... leaves 2 needs 1, 5 needs 1. Add (2,5).
    - $\langle 2 3 \rangle \langle 4 5 \rangle \langle 2 5 \rangle$.
    - (2,4), (3,5) ... leaves 2 needs 1, 5 needs 1. Add (2,5).
    - $\langle 2 4 \rangle \langle 3 5 \rangle \langle 2 5 \rangle$.
    - (2,3), (2,4) ... 3,4 satisfied. 2 needs 0. 5 needs 2. Impossible (must connect 5).
    - Wait, (3,5) and (4,5)?
    - $\langle 3 5 \rangle \langle 4 5 \rangle$. 5 has deg 2. 3,4 deg 1. 2 needs 2.
    - Multiply by $\langle 2 2 \rangle$? No, zero.
    - We need degree 2 in 2.
    - Must involve (2,X) twice.
    - If X=5, ok.
    - If X=3? (2,3). Then 3 satisfied. 4 needs 1. 5 needs 2. 2 needs 1.
      - Add (2,4). 4 satisfied. 2 sat. 5 needs 2.
      - Add (5,5)? No.
      - Must connect 5.
- So we can iterate these very few graphs manually.
- **This is the solution.**

## Plan
1.  Define basis for Weights $[0,0,2,1,1,2]$ on nodes $\{2,3,4,5\}$.
2.  Fit $N / (\langle 0 1 \rangle^8 D_{cyclic})$.
3.  The result should be simple integers.








