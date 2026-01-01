# CURRENT BCFW FORMULA (NEEDS CORRECTION)

## Current Implementation

**File**: `proper_proof_amplituhedron_hodges.sage`
**Function**: `amplituhedron_from_bcfw_cells()`

### Channels
For 6-point MHV gravity, we use 3+3 factorization channels:
- (0,3): particles 0,1,2 | 3,4,5
- (1,4): particles 1,2,3 | 4,5,0
- (2,5): particles 2,3,4 | 5,0,1
- (3,0): particles 3,4,5 | 0,1,2
- (4,1): particles 4,5,0 | 1,2,3
- (5,2): particles 5,0,1 | 2,3,4

### Current Formula (WRONG - PLACEHOLDER)

For each channel `(i, j)` where channel = particles `i, i+1, i+2`:

```python
# Channel 4-bracket
j0 = first outside particle
j0m1 = (j0 - 1) % n
channel_bracket = <i i+1 j0m1 j0>

# Channel angle brackets
channel_angles = <i i+1> * <i+1 i+2> * <i+2 i+3>

# All angle brackets (Parke-Taylor denominator)
all_angles = <0 1> * <1 2> * <2 3> * <3 4> * <4 5> * <5 0>

# CURRENT TERM (WRONG):
term = (channel_angles^2) / (channel_bracket^2 * all_angles^2)
```

### What We Know

1. **BCFW recursion**: M_6 = sum_{channels} M_3_left * M_3_right / s_channel
2. **For MHV**: M_3 = 1 (up to factors)
3. **Channel invariant s**: In momentum twistors, related to 4-brackets
4. **Gravity structure**: Needs squared Parke-Taylor factors

### What We Need

The correct formula that makes:
```
sum_{all channels} term = Hodges_determinant
```

### Issues with Current Formula

- Scale factor varies wildly (not constant)
- Doesn't match Hodges
- Powers/exponents may be wrong
- Missing factors or normalization

### Next Steps

1. Find correct BCFW formula from literature
2. Or derive from first principles
3. Or work backwards from Hodges to determine correct cell contributions

---

*Status: Placeholder formula - needs correction*









