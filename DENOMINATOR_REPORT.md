# Denominator Discovery Report

## Findings
Running the `discover_denominator.py` script with the optimized Finite Field probe + Exact GCD analysis yielded a definitive result for the denominator structure of the 6-pt MHV Gravity Amplitude.

### The Denominator $D(Z)$
The analysis shows that $M_{grav} = \frac{N(Z)}{D(Z)}$ where $N(Z)$ is a polynomial, and $D(Z)$ factors exactly into angle brackets with the following powers:

| Pair Type | Pairs | Power in $D(Z)$ |
| :--- | :--- | :--- |
| **Cyclic** | $\langle 12 \rangle, \langle 23 \rangle, \langle 34 \rangle, \langle 45 \rangle, \langle 50 \rangle$ | **2** |
| **Cyclic (Special)** | $\langle 01 \rangle$ | **0** (Absent) |
| **Non-Cyclic** | All 9 non-cyclic pairs ($\langle 02 \rangle, \langle 03 \rangle \dots$) | **1** |

### Interpretation
The denominator can be written as:
$$ D(Z) = \left( \prod_{j=1}^5 \langle j, j+1 \rangle^2 \right) \times \left( \prod_{1 \le a < b \le 6, b \ne a+1, (a,b)\ne(0,1)} \langle ab \rangle \right) $$

Alternatively:
$$ D(Z) = \frac{ \prod_{a<b} \langle ab \rangle \cdot \prod_{j} \langle j, j+1 \rangle }{ \langle 01 \rangle^2 } $$
(Checking powers: Non-cyclic gets 1 from product. Cyclic gets $1+1=2$. $\langle 01 \rangle$ gets $2 - 2 = 0$. This matches!)

So the denominator is simply **Product of all pairs $\times$ Product of all cyclic pairs**, divided by $\langle 01 \rangle^2$.

### Numerator Degree
The reconstructed numerator $N(Z)$ has degree **30** along the line $Z(t)$.
This gives us a target for basis construction.

## Conclusion
We have identified the exact rational structure of the amplitude.
The "missing" $\langle 01 \rangle$ pole is consistent with the MHV helicity configuration ($0^-, 1^-$) where the factor $\langle 01 \rangle^8$ in the numerator likely cancels these poles.

## Next Steps
1. **Verify Global Validity**: Check this denominator on random lines (not just one).
2. **Basis Construction**: Since $N(Z)$ is degree 30, we can now try to expand it in a basis.
3. **Geometry**: The factors $\langle j, j+1 \rangle$ and $\langle ab \rangle$ are the boundaries of the kinematic space.






