# Known Bugs and Fixes

## Bug 1: Sage Integer in random.seed()

**Symptom**:
```
TypeError: The only supported seed types are: None, int, float, str, bytes, and bytearray.
```

**Cause**: Sage preparser converts `42` to `Integer(42)`. Python's `random.seed()` rejects Sage integers.

**Fix**: Wrap all seed values with `int()`:
```python
# Bad
random.seed(42)

# Good
random.seed(int(42))
```

**Files affected**: Multiple files in `src/`, `tests/`, `tools/`. Search for `random.seed` to find them.

---

## Bug 2: Module import failures in Sage

**Symptom**:
```
ModuleNotFoundError: No module named 'src'
```

**Cause**: Sage's `load()` doesn't set up Python import paths.

**Fix**: Use `load()` instead of `import`:
```python
# Bad
from src.hodges import hodges_6pt_mhv_reduced

# Good
load("src/hodges.sage")
# Now hodges_6pt_mhv_reduced is in global namespace
```

---

## Bug 3: SyntaxWarning about escape sequences

**Symptom**:
```
SyntaxWarning: invalid escape sequence '\d'
```

**Cause**: Docstrings containing LaTeX-like patterns (e.g., `\det`, `\bar`).

**Fix**: Use raw strings for docstrings with backslashes:
```python
# Bad
"""Formula: \det(M)"""

# Good  
r"""Formula: \det(M)"""
```

---

## Bug 4: Empty intersection in DCP search

**Symptom**: After intersecting with first boundary, dimension drops to 0.

**Possible causes**:
1. S6 symmetry is too restrictive
2. Chart choice is incompatible with boundary
3. Boundary constraint matrix has numerical issues

**Debugging steps**:
1. Print the invariant basis vectors explicitly
2. Evaluate each basis vector on boundary kinematics
3. Check if any vector survives (should have rank < dim)

---

## Bug 5: KLT/Hodges ratio varies

**Symptom**: `M_KLT / M_Hodges` is not constant across samples.

**Explanation**: This is NOT a bug. The two formulas compute amplitudes with different little-group weights:
- Hodges Reduced: weight -4 per particle (total -24)
- KLT: weight +4 per particle (total +24 → -8 after PT factors)

**Resolution**: Use Hodges Reduced as the oracle. It's the "stripped" amplitude without helicity factors.

