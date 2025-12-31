# Physics Context for Non-Experts

## What is this project about?

We're trying to prove that **gravity scattering amplitudes** can be derived purely from **geometry** — without using Feynman diagrams or field theory.

## The Physical Object

**6-point MHV Gravity Amplitude**: The probability amplitude for 6 gravitons to scatter, in the "Maximal Helicity Violating" configuration (2 negative helicity, 4 positive helicity).

Two equivalent formulas exist:
1. **Hodges formula**: A determinant-based expression
2. **KLT formula**: "Gravity = (Yang-Mills)² " double-copy construction

## The Mathematical Goal

Show that the Hodges amplitude is the **unique** differential form satisfying:
1. **Pole structure**: Only poles at `s_S = 0` (where `s_S` is a Mandelstam invariant)
2. **Symmetry**: Invariant under permutations of external particles
3. **Factorization**: Residues at poles match products of lower-point amplitudes

If proven, this would mean: **consistency alone determines the amplitude**.

## Key Concepts

### Spinor Helicity Variables
Momenta are written as products of 2-component spinors:
```
p_i^μ = λ_i^α * λ̃_i^α̇
```
Angle brackets: `<ij> = det(λ_i, λ_j)`
Square brackets: `[ij] = det(λ̃_i, λ̃_j)`

### Mandelstam Invariants
```
s_ij = (p_i + p_j)² = <ij>[ij]
s_ijk = (p_i + p_j + p_k)² = s_ij + s_jk + s_ik
```

### Momentum Conservation
```
Σ p_i = 0  →  Σ λ_i λ̃_i = 0
```
This constrains the kinematics to a (3n-10)-dimensional space.

### Positive Geometry
A geometric region whose "canonical form" (a differential form with logarithmic singularities on boundaries) equals a physical amplitude. The amplituhedron is the famous example for Yang-Mills.

## Why This Matters

If we can derive gravity amplitudes from geometry:
1. It suggests a deeper mathematical structure underlying quantum gravity
2. It could lead to new computational methods
3. It connects scattering amplitudes to algebraic geometry and combinatorics

## References

- Hodges, "A simple formula for gravitational MHV amplitudes" (arXiv:1204.1930)
- Arkani-Hamed & Trnka, "The Amplituhedron" (arXiv:1312.2007)
- Cachazo, He, Yuan (CHY), "Scattering equations" (arXiv:1306.6575)

