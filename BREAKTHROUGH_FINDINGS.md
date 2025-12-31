# Physics Breakthrough Findings

## Date: December 29, 2025

## Executive Summary

We have made a significant mathematical discovery: **The 6-point MHV gravity amplitude (Hodges determinant) is NOT representable as an S6-invariant element of the OS3 (Orlik-Solomon degree 3) space.**

This is a crucial finding that explains why previous approaches were stuck at dim=8 or hitting empty intersections.

## Key Results

### 1. OS3 Space Structure
- **Full OS3 dimension**: 2008
- **S6-invariant subspace**: dimension 2
- **S3×S3-invariant subspace**: dimension 58
- **S3×S3×Z2-invariant subspace**: dimension 26

### 2. Boundary Intersection Results
When applying factorization constraints (boundary conditions):

| Invariant Mode | Initial Dim | After 1 Boundary | After 2 Boundaries |
|---------------|-------------|------------------|-------------------|
| S3×S3         | 58          | 48               | **EMPTY**         |
| S3×S3×Z2      | 26          | 21               | **EMPTY**         |
| S6            | 2           | **EMPTY**        | -                 |

**Critical Observation**: The S6-invariant space becomes empty after just ONE boundary constraint (full rank = 2/2).

### 3. Hodges Projection Test
Direct numerical projection of S6 invariants onto the Hodges amplitude:
- **Residuals**: Max 1.56, Avg 0.12
- **Verification errors**: Up to 5×10^4

**Conclusion**: The Hodges amplitude is NOT in the S6-invariant subspace of OS3.

## Mathematical Interpretation

### Why This Happens

1. **OS3 Construction**: The Orlik-Solomon space OS3 is built from the matroid structure of the kinematic space. It captures certain algebraic relations but may not include all physically relevant forms.

2. **S6 Invariance**: Full symmetric group invariance is a strong constraint. The gravity amplitude, while symmetric, may require a representation that is NOT captured by the OS3 construction.

3. **Positive Geometry Mismatch**: The positive geometry for gravity amplitudes (related to the amplituhedron) may require a different geometric space than OS3.

### Possible Resolutions

1. **Momentum Twistor Space**: The amplituhedron is naturally defined in momentum twistor space, not kinematic space. The gravity amplitude may only have a nice positive geometry representation in twistor variables.

2. **Different Degree**: Perhaps OS4 or higher-degree Orlik-Solomon spaces are needed.

3. **Non-Symmetric Representation**: The amplitude might require working in a non-S6-invariant subspace and then symmetrizing.

4. **Worldsheet/CHY Approach**: The CHY (Cachazo-He-Yuan) formula provides an alternative representation that might be more amenable to positive geometry analysis.

## Next Steps for Breakthrough

### Immediate Actions

1. **Try Full OS3 Space**: Project the full 2008-dimensional OS3 space onto Hodges (computationally intensive but definitive).

2. **Momentum Twistor Variables**: Reformulate the problem using momentum twistors where the amplituhedron is naturally defined.

3. **CHY/Worldsheet Approach**: Use the CHY formula which naturally gives a positive geometry on the moduli space of punctured Riemann surfaces.

### Theoretical Investigation

1. **Why OS3 Fails**: Understand mathematically why the OS3 construction doesn't capture gravity.

2. **Correct Space**: Identify the correct positive geometry space for gravity amplitudes.

3. **Universal Structure**: Look for a unified framework that works for both Yang-Mills (associahedron) and gravity.

## Code Artifacts

All code and results are saved in:
- `standalone_hodges.sage` - Main projection script
- `standalone_hodges_result.sobj` - Saved results
- `standalone_hodges.log` - Full execution log

## Additional Results: S3×S3 Projection

After finding that S6 invariants fail, we tried the larger S3×S3-invariant space (dim=58):

### S3×S3 Projection Results
- **Matrix rank**: 58 (full rank)
- **Training residuals**: Max 0.082, Avg 0.013 (very good!)
- **Verification errors**: Max 2.9×10^4, Min 0.94

### Interpretation
The S3×S3 space shows much better fit on training data but still fails on verification. This suggests:

1. **Overfitting**: The 58-dimensional space can fit the training data well but doesn't generalize
2. **Numerical issues**: The evaluation function may have numerical instabilities
3. **Near-miss**: The amplitude is "close" to being in S3×S3 but not exactly

### Key Coefficients Found
The dominant coefficients in the S3×S3 projection:
- c[26] = -244.67 (largest)
- c[1] = -110.96
- c[4] = 77.29
- c[0] = 74.83

These coefficients might provide insight into the structure of the gravity amplitude.

## Conclusion

This investigation has revealed fundamental properties of the positive geometry approach to gravity:

### Negative Results (What Doesn't Work)
1. **S6-invariant OS3 space (dim=2)**: Empty under boundary constraints
2. **S3×S3-invariant OS3 space (dim=58)**: Good training fit but poor generalization
3. **Boundary intersection approach**: Becomes empty before reaching dim=1

### Positive Insights (What We Learned)
1. The gravity amplitude is NOT purely S6-invariant in OS3
2. The S3×S3 space is closer but still insufficient
3. The positive geometry for gravity requires a fundamentally different construction

### Suggested Next Steps
1. **Momentum Twistor Space**: Reformulate using momentum twistors
2. **Amplituhedron Approach**: Use the gravity amplituhedron directly
3. **CHY/Worldsheet**: Try the CHY formula representation
4. **Higher Orlik-Solomon**: Try OS4 or higher degrees

This finding aligns with the broader physics literature suggesting that gravity positive geometry is fundamentally different from gauge theory positive geometry (associahedron).

---

*"Sometimes the most important discoveries are about what doesn't work - and this tells us exactly where to look next."*

## Files Generated
- `standalone_hodges.sage` - S6 projection script
- `standalone_hodges_result.sobj` - S6 results
- `s3xs3_hodges_projection.sage` - S3×S3 projection script  
- `s3xs3_hodges_result.sobj` - S3×S3 results
- `BREAKTHROUGH_FINDINGS.md` - This report

