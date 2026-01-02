# Amplituhedron Implementation Plan

## Current Status

We have identified that:
1. ✅ OS3/S6 approach fails for gravity
2. ✅ Momentum twistor framework is the correct approach
3. ✅ Amplituhedron cells can be identified (9 cells for 6-point MHV)
4. ⚠️ Need to verify cell sum equals Hodges amplitude

## Implementation Strategy

### Phase 1: Basic Amplituhedron Structure ✅ (In Progress)

**Goal**: Understand the cell structure and compute basic forms

**Components**:
- [x] Momentum twistor kinematics
- [x] Angle brackets <ij>
- [x] 4-brackets <ijkl>
- [x] Cell enumeration
- [ ] Cell canonical forms (in progress)

### Phase 2: BCFW Cell Decomposition

**Goal**: Compute amplitude as sum over BCFW cells

**Components**:
- [ ] Proper BCFW recursion implementation
- [ ] Cell-to-cell transition rules
- [ ] Factorization on boundaries
- [ ] Sum over all cells

### Phase 3: Verification

**Goal**: Verify amplituhedron sum equals Hodges formula

**Components**:
- [ ] Compare on multiple kinematic points
- [ ] Check numerical accuracy
- [ ] Verify factorization properties
- [ ] Test soft limits

### Phase 4: Optimization and Extension

**Goal**: Make it robust and extend to higher points

**Components**:
- [ ] Optimize cell enumeration
- [ ] Cache computations
- [ ] Extend to 7-point, 8-point
- [ ] Handle NMHV cases

## Key Formulas

### Momentum Twistor Relations

1. **Angle brackets**: <ij> = Z_i^1 Z_j^2 - Z_i^2 Z_j^1
2. **Square brackets**: [ij] = <i-1 i j-1 j> / (<i-1 i> <j-1 j>)
3. **4-brackets**: <ijkl> = det(Z_i, Z_j, Z_k, Z_l)

### Hodges Formula in Twistors

For 6-point MHV:
```
M_6 = det'(Phi) / (<12><23><34><45><56><61>)
```

where Phi_{ij} = [ij]/<ij> for i≠j, and Phi_{ii} has diagonal terms.

### BCFW Terms

Each BCFW term corresponds to a factorization channel:
```
M_n = sum_{channels} M_L * 1/P^2 * M_R
```

In momentum twistors, channels are labeled by 4-brackets.

## Known Issues and Solutions

### Issue 1: Singular Kinematics
**Problem**: Some kinematic points give zero brackets
**Solution**: Filter out singular points, use generic kinematics

### Issue 2: Cell Form Computation
**Problem**: Need correct formula for canonical forms
**Solution**: Use known BCFW formulas, verify against literature

### Issue 3: Sign Conventions
**Problem**: 4-brackets have sign ambiguities
**Solution**: Use consistent ordering, check against known results

## Testing Strategy

1. **Unit Tests**: Test individual components (brackets, cells)
2. **Integration Tests**: Test full amplitude computation
3. **Verification Tests**: Compare with Hodges formula
4. **Robustness Tests**: Test on many random kinematic points

## Success Criteria

1. ✅ Can compute momentum twistor brackets
2. ✅ Can identify amplituhedron cells
3. ⚠️ Can compute cell canonical forms (in progress)
4. ⬜ Cell sum equals Hodges amplitude (within numerical precision)
5. ⬜ Works on multiple kinematic points
6. ⬜ Handles singular cases gracefully

## Next Immediate Steps

1. **Fix cell canonical form computation**
   - Use correct BCFW formula
   - Verify signs and factors
   - Test on known examples

2. **Implement proper BCFW recursion**
   - Handle 3-point amplitudes correctly
   - Compute channel invariants
   - Sum over all channels

3. **Add verification framework**
   - Compare with Hodges on many points
   - Check relative errors
   - Identify systematic issues

## Resources

### Papers
- Arkani-Hamed, Trnka - "The Amplituhedron" (arXiv:1312.2007)
- Hodges - "Eliminating spurious poles" (arXiv:0905.1473)
- Arkani-Hamed et al. - "All-Loop Integrand" (arXiv:1008.2958)

### Code References
- Amplituhedron Mathematica package
- BCFW recursion implementations
- Momentum twistor tools

## Timeline

- **Week 1**: Complete Phase 1 (cell structure) ✅
- **Week 2**: Complete Phase 2 (BCFW decomposition) ⬜
- **Week 3**: Complete Phase 3 (verification) ⬜
- **Week 4**: Complete Phase 4 (optimization) ⬜

---

*"The amplituhedron is the geometry of scattering amplitudes."*










