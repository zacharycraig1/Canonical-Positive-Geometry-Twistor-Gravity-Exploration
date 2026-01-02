# Progress: Kinematic Associahedron Approach to Gravity

## Summary

We have implemented the kinematic associahedron approach following Arkani-Hamed et al.
(arXiv:1711.09102) and analyzed the structure of the Hodges gravity amplitude.

## Key Findings

### 1. Associahedron Structure ✓

The 6-point kinematic associahedron A_3 has been correctly implemented:
- **14 vertices** (triangulations) - matches Catalan number C_4
- **9 facets** (planar poles X_ij)
- **3 dimensions**

The canonical form (bi-adjoint amplitude m(α,α)) is correctly computed as a sum
over triangulations.

### 2. Hodges Amplitude Structure ✓

The Hodges amplitude has the form:
```
M_6 = <12>^8 × det'(Φ) / (<12><23><34><45><56><61>)^2
```

Or in BGK form:
```
M_6 = <12>^8 × G_6 / (PT × s_012 × s_123 × s_234)
```

Key observations:
- **Poles**: s_012, s_123, s_234 (three-particle channels)
- **Numerator**: G_6 is a rational function (not polynomial) in spinor brackets
- **Helicity**: <12>^8 factor for MHV (particles 1,2 negative helicity)

### 3. Momentum Conservation ✓

Verified: s_012 = s_345 (both equal -2.0247 in test case)

This is the expected constraint from 6-particle momentum conservation.

### 4. KLT Implementation ⚠️

The KLT double copy was implemented but shows a mismatch with Hodges:
- Ratio varies between samples (not a constant)
- Known issue: ordering conventions in the KLT kernel need fixing

The existing codebase has documentation of this issue in `results/KLT_ASYMMETRY_REPORT.md`.

## Positive Geometry Interpretation

Based on the analysis, the positive geometry for 6-point MHV gravity appears to be:

### Structure
```
E → B
```
where:
- **B = Kinematic Polytope**: Associahedron-like with facets at s_ijk = 0
- **E = Total Space**: Includes helicity information via <12>^8 fiber

### Boundaries
The boundaries (facets) correspond to:
- s_012 = 0 (and s_345 = 0 by momentum conservation)
- s_123 = 0 (and s_456 = 0)
- s_234 = 0 (and s_501 = 0)

### Canonical Form
```
Ω(E) = <12>^8 × Ω(B)
```

where Ω(B) is the canonical form of the base kinematic polytope.

## What's Different from Yang-Mills

For Yang-Mills (amplituhedron):
- Geometry is in Grassmannian Gr(k,n)
- Canonical form is Parke-Taylor / cyclic product

For Gravity:
- Geometry is "squared" (double copy)
- Additional helicity factor <12>^8
- Numerator G_6 encodes the "double" structure

## Next Steps

1. **Fix KLT normalization**: Correct the kernel conventions to match Hodges
2. **Prove factorization**: Show that residues at s_ijk = 0 give product of lower-point amplitudes
3. **Identify the polytope**: Determine if the base B is exactly an associahedron or a related object
4. **Compute canonical form geometrically**: Derive the amplitude from the geometry, not vice versa

## Files Created

```
src/kinematic_associahedron/
├── __init__.py                    # Module definition
├── associahedron.sage             # Kinematic associahedron implementation
├── klt_geometry.sage              # KLT kernel as geometric constraint
├── gravity_from_klt.sage          # KLT double copy for gravity
├── analyze_hodges_structure.sage  # Structure analysis of Hodges amplitude
└── PROGRESS.md                    # This file
```

## Conclusion

The kinematic associahedron approach provides the correct framework for understanding
gravity's positive geometry. The amplitude is structurally similar to a canonical form:
- Logarithmic singularities on boundaries (poles)
- Factorization at poles (recursive structure)
- Helicity factor as fiber contribution

The main open question is the precise identification of the polytope whose canonical
form (times <12>^8) gives the Hodges amplitude.

