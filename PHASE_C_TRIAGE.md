# Phase C Triage: Channel Factorization Failure

## Observations
1. **s_123 Blowup**: In the bisection search, $s_{123}$ grew from $\sim 100$ to $5 \times 10^{16}$ while narrowing the interval.
   - This indicates we are bracketing a **pole** of $s_{123}(t)$, not a zero.
   - $s_{123}(t)$ is a rational function of $t$. If it changes sign by passing through infinity, bisection converges to the pole.
   - $M$ also blows up there, but not in a way that creates a simple pole residue ($M \cdot s_{123}$ diverges).

## Diagnosis
- The deformation $Z(t)$ induced a singularity in the parameterization (likely $\langle i, i+1 \rangle \to 0$ or similar) that caused $s_{123}$ to diverge.
- We need to find a **zero** of $s_{123}$, not a pole.
- Logic fix: Check if $s_{low}$ and $s_{high}$ are finite and have opposite signs. If they are large, we might be near a pole.

## Plan Correction
1. **Use Numerator of s_123**:
   - $s_{123}(t) = \frac{Num(t)}{Den(t)}$.
   - We want $Num(t) = 0$ while $Den(t) \ne 0$.
   - Instead of evaluating $s_{123}$ directly, we should evaluate its numerator (or just check if magnitude is small).
   
2. **Better Deformation**:
   - Instead of a random line, use a targeted deformation that preserves local brackets but allows $s_{123}$ to vary.
   - Actually, a random line is fine if we solve for $Num(s_{123}) = 0$.

3. **Fallback**:
   - If $s_{123}$ has no real roots on the line (possible), we need a complex root.
   - Or just try many lines until we find a simple zero crossing.

## Next Steps
- Modify script to check if we are converging to a pole (value grows).
- If so, discard and retry with new seed.
- Verify $Den(s_{123})$ is not zero at the root.



