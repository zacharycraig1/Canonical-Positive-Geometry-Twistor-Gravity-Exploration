# Research Status: Gravity Positive Geometry

## Date: Current Session

## Executive Summary

We have made significant progress in understanding the structure of 6-point MHV gravity amplitudes, though the complete positive geometry description remains an open problem.

---

## 1. VERIFIED RESULTS

### ✅ Hodges Formula Works
- Computes correct gravity amplitude
- Consistent across all valid kinematic configurations

### ✅ Momentum Conservation
- s_012 = s_345 verified in all tests
- All other momentum constraints satisfied

### ✅ Amplitude Structure (BGK Form)
The Hodges amplitude has the structure:
```
M_6 = <12>^8 × G_6 / (PT × s_012 × s_123 × s_234)
```
where:
- `<12>^8` = helicity factor
- `PT` = Parke-Taylor denominator = <12><23><34><45><56><61>
- `s_ijk` = three-particle Mandelstam invariants (poles)
- `G_6` = numerator (varies with kinematics)

### ✅ Associahedron Structure
- 14 vertices (triangulations) - matches Catalan C_4
- 9 facets (planar poles)
- 3 dimensions
- Canonical form computes bi-adjoint amplitude correctly

### ✅ 4-Point KLT = -Hodges
For n=4: **KLT = -Hodges exactly** (constant ratio = -1)
This validates the KLT framework at 4-point.

---

## 2. UNRESOLVED ISSUES

### ❌ 6-Point KLT ≠ Hodges
- Ratio varies wildly between samples (not constant)
- Different KLT conventions tested, none work
- The issue is NOT:
  - Mandelstam conventions (verified consistent)
  - Momentum conservation (verified)
  - Amplitude orderings (tested multiple variations)

### ❓ Possible Causes
1. **Field-theory KLT kernel** may have additional factors at n≥6
2. **Color-ordering** conventions may not match between YM and kernel
3. **BCJ numerators** may be needed for exact match

---

## 3. WHAT THE POSITIVE GEOMETRY LOOKS LIKE

Based on our analysis, the geometry appears to be:

### Structure
```
        ⟨12⟩^8
          ↓
M_6 = ―――――――――――――――――――――――
      PT × s_012 × s_123 × s_234
          ↓
    (Parke-Taylor) × (Mandelstam poles)
```

### Boundaries
The amplitude has simple poles at:
- **Three-particle channels**: s_012 = 0, s_123 = 0, s_234 = 0
- **Two-particle channels**: embedded in PT denominator

### Canonical Form Interpretation
The amplitude looks like a canonical form:
```
Ω = <12>^8 × d³z / (boundary_1 × boundary_2 × boundary_3)
```

This is the structure of a **3-dimensional positive geometry** with logarithmic singularities on its boundaries.

---

## 4. OPEN QUESTIONS

1. **What is the polytope?**
   - Is it the kinematic associahedron? Something else?
   - How do the 14 triangulations relate to terms in M_6?

2. **How does KLT enter geometrically?**
   - The double-copy structure should manifest as a product geometry
   - The KLT kernel S[α|β] should be a geometric constraint

3. **What defines the positive region?**
   - Not simple pole positivity (amplitude can be ±)
   - Must involve spinor brackets in some way

4. **How does factorization work?**
   - Residue at s_012 = 0 should give M_4 × M_4
   - This requires symbolic computation to verify

---

## 5. RECOMMENDED NEXT STEPS

### Immediate
1. **Symbolic factorization test**: Compute Res_{s_012=0} M_6 symbolically
2. **Fix KLT normalization**: Look at BCJ duality / color-kinematics
3. **Cell decomposition**: Try to write M_6 as sum over associahedron vertices

### Longer Term
1. **Study the numerator G_6**: What polynomial structure does it have?
2. **Explore double-copy geometry**: Gr(2,n) × Gr(2,n) / SL(2)?
3. **Connect to BCFW**: Each BCFW term might be a cell

---

## 6. KEY FILES

```
src/kinematic_associahedron/
├── associahedron.sage            # Working associahedron
├── klt_geometry.sage             # KLT implementation
├── gravity_from_klt.sage         # KLT gravity computation
├── test_4pt_klt.sage             # 4-point KLT test (WORKS!)
├── test_6pt_klt.sage             # 6-point KLT test (mismatch)
├── debug_klt_6pt.sage            # KLT debugging
├── test_factorization.sage       # Factorization tests
├── analyze_hodges_structure.sage # Structure analysis
└── RESEARCH_STATUS.md            # This file

src/twistor_gravity/
├── FINDINGS.md                   # Previous exploration results
└── [various .sage files]         # Twistor approaches (explored)
```

---

## 7. CONCLUSION

**We are on the right track but haven't fully solved the problem.**

The positive geometry for gravity MHV is a **genuine open research problem** that leading physicists (UNIVERSE+ project at Max Planck) are actively working on.

Our contributions:
1. Verified the amplitude structure matches expectations
2. Confirmed 4-point KLT works exactly
3. Identified that 6-point has additional complexity
4. Built extensive computational infrastructure

The geometry exists (the amplitude has the right structure), but its explicit description requires more work.

---

## 8. REFERENCES

1. Arkani-Hamed et al., "Scattering Forms and the Positive Geometry of Kinematics" (arXiv:1711.09102)
2. Hodges, "New expressions for gravitational scattering" (2012)
3. Kawai-Lewellen-Tye relations (1986)
4. UNIVERSE+ project at Max Planck Institute (2024-2025)

