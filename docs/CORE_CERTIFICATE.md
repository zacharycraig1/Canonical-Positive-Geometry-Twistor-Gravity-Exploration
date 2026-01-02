# Core Certificate of Verification (Phase U)

This document certifies the reproducibility of the core identity:
**MHV Gravity Amplitude = Rooted Forest Polynomial evaluated on Spinor Edge Variables**

## Verified Statements
1. For n=4, 5, 6, the reconstructed MHV amplitude (from Hodges determinant / KLT) matches EXACTLY the forest polynomial evaluation times the known prefactors.
2. The prefactors are:
   \[
   M_{MHV} = (-1)^{n-1} \frac{\langle 01 \rangle^8}{(\prod_{k \notin R} C_k^2) (\prod_{cyc R} \langle r_i r_{i+1} \rangle^2)} F_{n,R}(z_{ij})
   \]
   where $z_{ij}$ are the spinor edge variables defined in `src/posgeom/physics_map.py`.

## Reproducibility Instructions
To verify the core identity, run the consolidated test suite from the project root:

```bash
python run_core_tests.py
```

This script will execute:
- `src/scripts/physics_pullback_n4.sage`
- `src/scripts/physics_pullback_n5.sage`
- `src/scripts/physics_pullback_n6.sage`

### Expected Output
The script should output "SUCCESS" for all checks.

### Environment
- Requires SageMath environment.
- Verified on Windows/Linux with SageMath 9.x/10.x.




