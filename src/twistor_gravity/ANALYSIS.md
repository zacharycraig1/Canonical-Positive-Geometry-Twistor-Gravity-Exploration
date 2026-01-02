# Twistor Space Gravity: Analysis of Approaches

## Executive Summary

This document analyzes three approaches to finding the positive geometry for 6-point MHV gravity amplitudes:

1. **Momentum Twistors / Gravituhedron**: Extend amplituhedron infrastructure
2. **Twistor String**: Skinner's worldsheet-to-twistor formula  
3. **Celestial Holography**: Map to 2D CFT on celestial sphere

Each approach offers different perspectives on positivity and the canonical form structure.

---

## 1. Current Status: CHY Approach

### What We Found

The worldsheet (CHY) approach was our original path:
- **Formula**: `M_n = <12>^8 × Σ [Pf'(Ψ)]² / det'(Φ)`
- **Positive region attempt**: `z_4 > z_5 > z_6 > 0`

### Issues Discovered

1. **Normalization discrepancy**: Factor of ~2220x between CHY sum and Hodges
2. **Sign mismatch**: CHY gives positive, Hodges gives negative for same kinematics
3. **No simple positive region**: The 6 scattering equation solutions don't all lie in a single positive ordering chamber

### Key Insight

The CHY formula works but the "positive geometry" isn't visible in worldsheet moduli space M_{0,6}. The geometry may live in a different space.

---

## 2. Approach A: Momentum Twistors / Gravituhedron

### Theory

**Momentum twistors** Z_i ∈ CP^3 encode kinematics with automatic momentum conservation.

**Positive Grassmannian** Gr_+(k, n): All ordered minors positive.

For MHV (k=2), the amplituhedron A_{n,2} ⊂ Gr(2, n) gives Yang-Mills amplitudes.

### Hypothesis for Gravity

> Gravity geometry = (Amplituhedron) × (Amplituhedron) / SL(2)

This "double copy" structure mirrors the KLT relations:
```
M_gravity = Σ S[α|β] A_YM[α] A_YM[β]
```

### Positivity Condition

All ordered 4-brackets must be positive:
```
⟨i i+1 j j+1⟩ > 0  for all cyclically ordered i < i+1 < j < j+1
```

### Implementation

- `hodges_twistor.sage`: Hodges formula in pure twistor variables
- `gravituhedron.sage`: Candidate geometry with double-copy structure

### Pros
- Natural positivity conditions on Grassmannian
- Momentum conservation automatic
- Connects to known amplituhedron technology

### Cons
- Double-copy geometry not yet fully understood
- May need new mathematical structures beyond Grassmannian

---

## 3. Approach B: Twistor String

### Theory

**Key insight**: For MHV, all external twistors lie on a **line** in CP^3.

The worldsheet Σ (genus 0, n punctures) maps to this twistor line with degree 1.

### Formula Structure

```
M_n = ∫_{M_{0,n}} δ(degree-1 constraint) × [integrand]
```

The delta functions impose that the map is degree 1 (MHV condition).

### Relation to CHY

CHY localizes the integral to (n-3)! points via scattering equations.
Twistor string provides the same localization but with geometric meaning.

### Positivity

Positivity might manifest as:
- The twistor line lies in positive region of PT
- The worldsheet map has positive moduli

### Implementation

- `twistor_string.sage`: Line geometry and MHV vertex structure

### Pros
- Geometric meaning of MHV constraint
- Natural connection to twistor theory literature
- Proven equivalence to CHY

### Cons
- Full gravity formula requires careful treatment of integrand
- Positivity of "degree 1 map" not fully characterized

---

## 4. Approach C: Celestial Holography

### Theory

**Celestial amplitudes**: Mellin transform of momentum-space amplitudes.

```
Ã(Δ_i, z_i) = ∫ dω ω^{Δ-1} M(ω q_i)
```

The result is a 2D CFT correlator on the celestial sphere (sky at null infinity).

### Symmetry

MHV gravity amplitudes have **Lw_{1+∞}** symmetry (extension of BMS group).

Soft graviton theorems ↔ Ward identities for this algebra.

### Positivity

Potential positivity sources:
1. **OPE coefficients** may have definite sign
2. **Conformal block expansion** may be positive for unitary CFT
3. **Unitarity bounds**: Δ ≥ 1 for massless particles

### Implementation

- `celestial_map.sage`: Celestial coordinates and correlator structure

### Pros
- Fresh perspective on amplitude structure
- Deep connection to symmetries (soft theorems)
- May reveal hidden positivity via CFT techniques

### Cons
- Most speculative approach
- Full theory still under development
- Connection to momentum-space positivity unclear

---

## 5. Comparison Table

| Aspect | Momentum Twistor | Twistor String | Celestial |
|--------|------------------|----------------|-----------|
| **Space** | Gr_+(4, n) | Line in PT | Celestial S² |
| **Positivity** | ⟨i i+1 j j+1⟩ > 0 | Degree-1 positivity | OPE signs |
| **Formula** | Hodges det | MHV vertex² | CFT correlator |
| **Maturity** | Well-developed | Proven | Emerging |
| **Geometry** | Grassmannian | Moduli of lines | Unknown |
| **Double copy** | Built-in | From vertex | Via KLT |

---

## 6. Recommended Path Forward

### Phase 1: Validate Cross-Equivalence

Run `compare_all.sage` to verify all methods give the same amplitude:
```bash
sage compare_all.sage
```

Expected: Hodges (spinor) = Hodges (twistor) = Gravituhedron = Twistor string

### Phase 2: Focus on Momentum Twistors

This approach is most developed and has clearest positivity conditions.

**Key questions**:
1. What is the "gravituhedron" as a mathematical object?
2. How does the product/fiber structure work?
3. What are the boundaries (factorization channels)?

### Phase 3: Twistor String for Geometric Insight

Use twistor string to understand **why** positivity fails on worldsheet:
- The geometry is in twistor space, not moduli space
- MHV constraint (degree 1) is the key

### Phase 4: Celestial for Future Directions

Keep celestial as a long-term research direction:
- May provide new positivity constraints via CFT
- Connection to asymptotic symmetries and holography

---

## 7. Files Created

```
src/twistor_gravity/
├── __init__.py           # Module documentation
├── hodges_twistor.sage   # Hodges in twistor variables
├── gravituhedron.sage    # Candidate positive geometry
├── twistor_string.sage   # Skinner formula skeleton
├── celestial_map.sage    # Celestial amplitude transform
├── compare_all.sage      # Cross-validation framework
└── ANALYSIS.md           # This document
```

---

## 8. Key References

1. **Amplituhedron**: Arkani-Hamed & Trnka, arXiv:1312.2007
2. **Twistor String Gravity**: Skinner, arXiv:1301.0868
3. **Celestial Holography**: Pasterski et al., various 2024 papers
4. **Hodges Formula**: Hodges, arXiv:1204.1930
5. **CHY Formula**: Cachazo-He-Yuan, arXiv:1307.2199

---

## 9. Conclusion

The search for gravity's positive geometry is open. Three promising paths exist:

1. **Momentum twistors**: Most concrete, best tools
2. **Twistor string**: Geometric insight, proven equivalence
3. **Celestial**: Novel, may reveal hidden structure

The key insight from our CHY investigation: **the geometry isn't in worldsheet moduli space**. It likely lives in twistor space or a product thereof.

Next step: Run cross-validation to establish baseline, then focus on characterizing the gravituhedron boundaries.

