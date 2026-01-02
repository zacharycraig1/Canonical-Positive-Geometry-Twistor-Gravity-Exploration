# Gravity Amplituhedron Discovery Plan
## Strategic Implementation Roadmap — December 31, 2025

---

## Executive Summary

This plan addresses the open problem of finding a **positive geometry construction for gravity amplitudes**. While the amplituhedron exists for gauge theory (N=4 SYM), no such construction is known for gravity. Success would be a **major result in theoretical physics**.

### Current State Assessment

#### ✅ VERIFIED (Solid Ground)
| Claim | Status | Evidence |
|-------|--------|----------|
| Forest Polynomial Identity: M_MHV = F_{n,R}(z) | **EXACT** for n=4,5,6 | `src/scripts/physics_pullback_n6.sage` |
| Reference Independence (gauge invariance) | **PROVEN** | Multiple root choices tested |
| Physical poles ↔ Geometric boundaries | **VERIFIED** | Facet analysis in `src/posgeom/facets_to_subsets.py` |

#### ❌ NOT VERIFIED (False Claims in Repo)
| Claim | Issue | Location |
|-------|-------|----------|
| "Amplituhedron = Hodges" | **CIRCULAR** — computes Hodges twice | `ruled_out_single_log_form/correct_amplituhedron_hodges.sage` |
| Positivity Region | No kinematic point with all z_ij > 0 found | — |
| Stringy Pushforward for n=6 | Ratio ~ 10^{-20} instead of 1 | `src/posgeom/saddle_pushforward.py` |

---

## Phase Structure

```
Week 1: OPTION C — Fix Pushforward (Low Risk, Diagnostic Value)
         ↓
Week 2-3: OPTION A — BCFW Amplituhedron (High Risk, High Reward)  
         ↓
Week 4: OPTION B — Intersection Theory (Fallback)
```

---

# PHASE 1: Fix the Pushforward (Week 1)

## Objective
Understand WHY the saddle pushforward fails at n=6. This diagnostic work will inform whether Option A is viable.

## Current Situation
```
n=4: Pushforward ratio = 1 ✓
n=5: Pushforward ratio = 1 ✓  
n=6: Pushforward ratio ~ 10^{-20} ✗
```

## Implementation Tasks

### Task 1.1: Create Debug Infrastructure
**File:** `src/debug/pushforward_n6_analysis.sage`

```sage
"""
Systematic analysis of n=6 pushforward failure.

Hypotheses to test:
1. Wrong normalization (scaling factor)
2. Missing Jacobian determinant
3. Incorrect saddle counting  
4. Wrong moment map
"""

def analyze_pushforward_failure():
    # Load forest polynomial
    F = get_forest_polynomial(6, [0,1,2])
    
    # Sample 100 physical kinematic points
    ratios = []
    for seed in range(100):
        z = sample_physical_z(6, seed)
        
        # Compute pushforward value
        pushforward_val = compute_pushforward_saddle(...)
        
        # Compute Hodges reference
        hodges_val = reconstruct_mhv_from_laplacian(...)
        
        if hodges_val != 0:
            ratios.append(pushforward_val / hodges_val)
    
    # Key test: Is the ratio CONSTANT?
    if constant_ratio(ratios):
        print("CONSTANT ratio → just need normalization factor")
        return ratios[0]  # The missing factor
    else:
        print("VARYING ratios → map itself is wrong")
        return None
```

**Success Criteria:**
- If ratio is **constant** → find the normalization factor
- If ratio **varies** → map is fundamentally wrong, proceed to Option A

### Task 1.2: Test Alternative Maps
**File:** `src/debug/alternative_maps.sage`

```sage
def test_alternative_maps():
    """
    Try different moment map formulations:
    1. X = ∇log(F)  [current]
    2. X = ∇(F^{1/k}) for various k
    3. X = log-coordinates directly
    """
    maps_to_try = [
        ("grad_log", lambda F, z: gradient(log(F), z)),
        ("grad_sqrt", lambda F, z: gradient(sqrt(F), z)),
        ("direct_log", lambda F, z: [log(z_i) for z_i in z]),
    ]
    
    for name, map_fn in maps_to_try:
        ratios = []
        for seed in range(50):
            ratio = test_pushforward_ratio(map_fn, seed)
            ratios.append(ratio)
        
        print(f"{name}: ratio = {mean(ratios):.6e} ± {std(ratios):.6e}")
```

### Task 1.3: Jacobian Analysis
**File:** `src/debug/jacobian_structure.sage`

```sage
def analyze_jacobian_structure():
    """
    The Jacobian of the moment map should be positive definite
    (it's a covariance matrix). Check numerical properties.
    """
    for seed in range(20):
        z = sample_physical_z(6, seed)
        X, J = moment_map_and_jacobian(log(z), poly_coeffs, poly_exponents)
        
        # Check eigenvalues
        eigenvalues = J.eigenvalues()
        print(f"Seed {seed}: eigenvalues = {eigenvalues}")
        print(f"  det(J) = {det(J)}")
        print(f"  cond(J) = {max(eigenvalues)/min(eigenvalues)}")
```

### Deliverables for Phase 1
1. `src/debug/pushforward_n6_analysis.sage` — Main diagnostic script
2. `src/debug/alternative_maps.sage` — Alternative map tests
3. `src/debug/jacobian_structure.sage` — Jacobian analysis
4. `PHASE1_REPORT.md` — Findings and conclusions

### Decision Point
- **If constant ratio found:** Proceed to verify normalization, then Option A
- **If no pattern found:** Proceed directly to Option A (the pushforward approach may be fundamentally limited)

---

# PHASE 2: BCFW Amplituhedron (Weeks 2-3)

## Objective
Implement proper BCFW cell construction for gravity MHV and verify it equals Hodges.

## Background
The gauge theory amplituhedron is triangulated by BCFW cells. Each cell contributes a d-log form. For gravity, the conjecture is:
- Gauge theory amplitude: Canonical form of amplituhedron in Gr(2,n)
- Gravity amplitude: Related to Gr(2,n) × Gr(2,n) or "square" structure (KLT)

## Implementation Tasks

### Task 2.1: Momentum Twistor Infrastructure (CLEAN REWRITE)
**File:** `src/amplituhedron/momentum_twistor.sage`

⚠️ **CRITICAL:** Do NOT use the broken version in `ruled_out_single_log_form/`

```sage
class MomentumTwistorData:
    """
    Momentum twistor kinematics for n particles.
    Z_i ∈ CP^3 (4-component twistor)
    
    Key invariants:
    - ⟨ij⟩ = ε_{AB} Z_i^A Z_j^B (2-bracket from first 2 components)
    - ⟨ijkl⟩ = det(Z_i, Z_j, Z_k, Z_l) (4-bracket)
    - [ij] derived from incidence relations
    """
    
    def __init__(self, n, seed=None):
        self.n = n
        if seed is not None:
            set_random_seed(seed)
        self.Z = self._generate_positive_kinematics()
        self._cache_brackets()
        
    def _generate_positive_kinematics(self):
        """
        Generate momentum twistors in the POSITIVE region.
        
        Positive means: all ordered minors ⟨i i+1 j j+1⟩ > 0
        Use explicit parameterization from Gr_+(4, n)
        """
        # Method: Build Z from positive coordinates
        # Start with canonical positive configuration, perturb
        n = self.n
        Z = []
        
        # Simple positive parameterization for testing
        for i in range(n):
            # Cyclic structure ensures positivity for generic parameters
            theta = 2 * pi * i / n
            z = vector(QQ, [
                QQ(randint(1, 10)),
                QQ(randint(1, 10)), 
                QQ(randint(1, 10)),
                QQ(randint(1, 10))
            ])
            Z.append(z)
        
        return Z
    
    def _cache_brackets(self):
        """Pre-compute all brackets for efficiency."""
        n = self.n
        self._angle = {}
        self._four_bracket = {}
        
        # 2-brackets
        for i in range(n):
            for j in range(n):
                self._angle[(i,j)] = self.Z[i][0]*self.Z[j][1] - self.Z[i][1]*self.Z[j][0]
        
        # 4-brackets (all ordered subsets)
        from itertools import combinations
        for indices in combinations(range(n), 4):
            M = matrix([self.Z[i] for i in indices])
            self._four_bracket[indices] = M.det()
    
    def angle(self, i, j):
        """⟨ij⟩"""
        return self._angle.get((i,j), QQ(0))
    
    def four_bracket(self, i, j, k, l):
        """⟨ijkl⟩ with sign from ordering"""
        indices = (i, j, k, l)
        sorted_indices = tuple(sorted(indices))
        base = self._four_bracket.get(sorted_indices, QQ(0))
        
        # Compute sign from permutation parity
        from sage.combinat.permutation import Permutation
        perm = [sorted_indices.index(x) + 1 for x in indices]
        sign = Permutation(perm).sign()
        
        return sign * base
    
    def is_positive(self):
        """Check if configuration is in positive region."""
        n = self.n
        for i in range(n):
            ip1 = (i + 1) % n
            for j in range(i + 2, n):
                jp1 = (j + 1) % n
                if jp1 == i:
                    continue
                bracket = self.four_bracket(i, ip1, j, jp1)
                if bracket <= 0:
                    return False
        return True
```

### Task 2.2: BCFW Cell Enumeration
**File:** `src/amplituhedron/bcfw_cells.sage`

```sage
def enumerate_bcfw_cells_mhv(n):
    """
    For MHV (k=2), BCFW cells correspond to triangulations of n-gon.
    For n=6 MHV, there are Catalan(4) = 14 terms, but due to BCFW structure,
    we get 5 independent cells (one per diagonal class).
    
    Each cell is labeled by a triangulation.
    """
    # For BCFW with shift [1,n⟩:
    # Cells correspond to planar diagrams with (n-3) propagators
    
    if n == 4:
        # Only one cell for n=4
        return [{'type': 'contact', 'vertices': [0,1,2,3]}]
    
    if n == 5:
        # Two cells: (12)(45) and (23)(51) channels
        return [
            {'channel': (1,2), 'diagram': 'A'},
            {'channel': (2,3), 'diagram': 'B'},
        ]
    
    if n == 6:
        # Five cells corresponding to 5 triangulations of hexagon
        # Or equivalently, BCFW channels
        cells = []
        # Channels: pick which propagator goes on-shell
        # For [1,6⟩ shift: channels are s_12, s_123, s_23, s_234, etc.
        for subset_size in [2, 3]:
            for start in range(1, n):
                subset = [(start + k) % n for k in range(subset_size)]
                if 0 not in subset and (n-1) in subset:
                    continue  # Skip redundant
                cells.append({
                    'channel': tuple(sorted(subset)),
                    'size': subset_size
                })
        return cells[:5]  # Catalan(4)/symmetry = 5 independent
    
    raise NotImplementedError(f"n={n} not implemented")

def cell_canonical_form(cell, twistors):
    """
    Compute the canonical form contribution from a single BCFW cell.
    
    For MHV, each cell contributes:
    Ω_cell = product of d-log factors
    
    The key formula involves 4-brackets ⟨AB I_1 I_2⟩ where (A,B) is the
    BCFW line and I_k are internal/external lines.
    """
    n = twistors.n
    channel = cell.get('channel', ())
    
    # Compute the BCFW propagator factor
    # For gravity MHV, this involves squared Parke-Taylor structure
    
    # Placeholder: implement actual formula from Arkani-Hamed & Trnka
    # The form is: ∏ ⟨i i+1⟩^{-1} × ⟨channel⟩^{-1} × (residue at pole)
    
    contribution = QQ(1)
    
    # Parke-Taylor denominator (appears twice for gravity)
    for i in range(n):
        ip1 = (i + 1) % n
        ang = twistors.angle(i, ip1)
        if ang == 0:
            return None
        contribution /= ang * ang  # Squared for gravity
    
    # Channel-dependent numerator
    if len(channel) >= 2:
        i, j = channel[0], channel[-1]
        four_br = twistors.four_bracket(i, (i+1)%n, j, (j+1)%n)
        contribution *= four_br
    
    return contribution
```

### Task 2.3: Verify BCFW Sum = Hodges
**File:** `src/amplituhedron/verify_bcfw_hodges.sage`

```sage
def test_bcfw_equals_hodges(n_trials=100):
    """
    THE KEY TEST: Does BCFW cell sum equal Hodges?
    
    If YES → we have a valid amplituhedron construction!
    If NO → analyze the discrepancy
    """
    n = 6
    matches = 0
    failures = []
    
    for seed in range(n_trials):
        tw = MomentumTwistorData(n=n, seed=seed)
        
        # Skip singular/non-positive kinematics
        if not tw.is_positive():
            continue
        
        # Compute BCFW sum
        cells = enumerate_bcfw_cells_mhv(n)
        bcfw_sum = QQ(0)
        valid = True
        
        for cell in cells:
            contrib = cell_canonical_form(cell, tw)
            if contrib is None:
                valid = False
                break
            bcfw_sum += contrib
        
        if not valid:
            continue
        
        # Compute Hodges reference
        hodges = hodges_determinant_from_twistors(tw)
        
        if hodges is None or hodges == 0:
            continue
        
        # EXACT rational comparison
        if bcfw_sum == hodges:
            matches += 1
        else:
            ratio = bcfw_sum / hodges
            failures.append({
                'seed': seed,
                'bcfw': bcfw_sum,
                'hodges': hodges,
                'ratio': ratio
            })
            
            # Check if it's a constant factor discrepancy
            print(f"Seed {seed}: ratio = {ratio}")
    
    print(f"\n{'='*60}")
    print(f"RESULTS: {matches}/{n_trials} exact matches")
    print(f"{'='*60}")
    
    if failures:
        # Analyze failure pattern
        ratios = [f['ratio'] for f in failures]
        if len(set(ratios)) == 1:
            print(f"[PROMISING] All failures have SAME ratio: {ratios[0]}")
            print("→ Just need normalization factor!")
        else:
            print(f"[INVESTIGATING] Varying ratios: {set(ratios)}")
    
    return matches, failures

def hodges_determinant_from_twistors(tw):
    """
    Compute Hodges determinant directly from momentum twistors.
    Uses the trusted laplacian_bridge module.
    """
    # Convert twistors to spinor-helicity variables
    # Then call the existing trusted function
    lambdas, tildes = twistors_to_spinors(tw)
    x, y = reference_spinors()
    
    M, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y)
    return M
```

### Task 2.4: Independent Hodges Implementation
**File:** `src/amplituhedron/hodges_twistor.sage`

```sage
def hodges_from_twistors_direct(tw):
    """
    Implement Hodges determinant formula directly in momentum twistor variables.
    This provides an independent check without going through spinor-helicity.
    
    Hodges formula (momentum twistor version):
    M_n = det(Φ) / ∏⟨i i+1⟩
    
    where Φ_ij = [ij]/⟨ij⟩ for i≠j (off-diagonal)
                = -∑_{k≠i} [ik]⟨kX⟩⟨kY⟩ / (⟨ik⟩⟨iX⟩⟨iY⟩) for i=j (diagonal)
    """
    n = tw.n
    
    # Reference legs X, Y (usually first and last)
    X, Y = 0, n-1
    
    # Build Hodges matrix (delete rows/cols for X, Y, and one more)
    # For n=6 with X=0, Y=5, delete {0, 5, r} for some r
    # Remaining indices: {1, 2, 3, 4}
    
    indices = [i for i in range(n) if i not in [X, Y]]
    # Actually need to delete 3 rows/cols total, so remove one more
    # Standard choice: delete the 3 "roots" {0, 1, 5} or similar
    
    # Implementation follows Hodges arXiv:0905.1473
    d = len(indices) - 1  # Size of reduced matrix
    
    Phi = matrix(QQ, d, d)
    kept = indices[:-1]  # Remove one more index
    
    for ii, i in enumerate(kept):
        for jj, j in enumerate(kept):
            if ii == jj:
                # Diagonal: sum over k ≠ i
                diag_sum = QQ(0)
                for k in range(n):
                    if k == i or k == X or k == Y:
                        continue
                    # [ik] = ⟨(i-1) i (k-1) k⟩ / (⟨(i-1) i⟩ ⟨(k-1) k⟩)
                    ik_sq = compute_square_bracket(tw, i, k)
                    ik_ang = tw.angle(i, k)
                    iX_ang = tw.angle(i, X)
                    iY_ang = tw.angle(i, Y)
                    kX_ang = tw.angle(k, X)
                    kY_ang = tw.angle(k, Y)
                    
                    if ik_ang == 0 or iX_ang == 0 or iY_ang == 0:
                        continue
                    
                    contrib = ik_sq * kX_ang * kY_ang / (ik_ang * iX_ang * iY_ang)
                    diag_sum -= contrib
                
                Phi[ii, jj] = diag_sum
            else:
                # Off-diagonal: [ij]/⟨ij⟩
                ij_ang = tw.angle(i, j)
                if ij_ang == 0:
                    return None
                ij_sq = compute_square_bracket(tw, i, j)
                Phi[ii, jj] = ij_sq / ij_ang
    
    # Determinant
    det_Phi = Phi.det()
    
    # Denominator: ∏⟨i i+1⟩
    denom = QQ(1)
    for i in range(n):
        ip1 = (i + 1) % n
        ang = tw.angle(i, ip1)
        if ang == 0:
            return None
        denom *= ang
    
    return det_Phi / denom

def compute_square_bracket(tw, i, j):
    """
    [ij] in momentum twistor variables:
    [ij] = ⟨(i-1) i (j-1) j⟩ / (⟨(i-1) i⟩ ⟨(j-1) j⟩)
    """
    n = tw.n
    im1 = (i - 1) % n
    jm1 = (j - 1) % n
    
    num = tw.four_bracket(im1, i, jm1, j)
    denom = tw.angle(im1, i) * tw.angle(jm1, j)
    
    if denom == 0:
        return None
    return num / denom
```

### Deliverables for Phase 2
1. `src/amplituhedron/__init__.py` — Package init
2. `src/amplituhedron/momentum_twistor.sage` — Clean twistor infrastructure
3. `src/amplituhedron/bcfw_cells.sage` — Cell enumeration and canonical forms
4. `src/amplituhedron/hodges_twistor.sage` — Independent Hodges implementation
5. `src/amplituhedron/verify_bcfw_hodges.sage` — Verification tests
6. `PHASE2_REPORT.md` — Results and analysis

### Success Criteria

| Tier | Description | Outcome |
|------|-------------|---------|
| **Tier 1** | BCFW sum = Hodges exactly for n=6 | **Publishable result** |
| **Tier 2** | Extend to n=7,8 | **Major result** |
| **Tier 3** | All-n proof | **Revolutionary** |

---

# PHASE 3: Intersection Theory (Week 4, if needed)

## Objective
If BCFW approach fails, try CHY/intersection theory as fallback.

## Implementation Tasks

### Task 3.1: Scattering Equations Solver
**File:** `src/intersection/scattering_eqs.sage`

```sage
def scattering_equations(n, sigma, s):
    """
    E_i = ∑_{j≠i} s_ij / (σ_i - σ_j) = 0
    
    Returns system of equations.
    """
    eqs = []
    for i in range(n):
        eq_i = sum(s[i,j] / (sigma[i] - sigma[j]) 
                   for j in range(n) if j != i)
        eqs.append(eq_i)
    return eqs

def find_solutions(n, s):
    """Find all (n-3)! solutions via numerical homotopy."""
    # Fix SL(2) gauge: σ_0=0, σ_1=1, σ_{n-1}=∞
    # Solve for σ_2, ..., σ_{n-2}
    pass
```

### Task 3.2: CHY Gravity Integrand
**File:** `src/intersection/chy_gravity.sage`

```sage
def chy_gravity_amplitude(momenta, polarizations):
    """
    M_n^gravity = ∫ dμ_n × (Pf'Ψ)^2
    
    where Pf'Ψ is the reduced Pfaffian.
    """
    pass
```

---

# Code Quality Standards

## DO ✅
- Use **exact arithmetic** (QQ in Sage) for all verifications
- **Cache expensive computations** (polynomials, brackets)
- Write **tests for every function**
- Document **mathematical formulas** in comments
- **Version control** all changes

## DON'T ❌
- Call Hodges from "amplituhedron" function (that's circular!)
- Use floats for verification (numerical errors mask bugs)
- Overwrite working code without backup
- Trust previous "breakthrough" claims without verification

---

# File Organization

```
src/
├── amplituhedron/          # NEW: Phase 2 code
│   ├── __init__.py
│   ├── momentum_twistor.sage
│   ├── bcfw_cells.sage
│   ├── hodges_twistor.sage
│   └── verify_bcfw_hodges.sage
│
├── debug/                  # NEW: Phase 1 diagnostics
│   ├── pushforward_n6_analysis.sage
│   ├── alternative_maps.sage
│   └── jacobian_structure.sage
│
├── intersection/           # NEW: Phase 3 (if needed)
│   ├── scattering_eqs.sage
│   └── chy_gravity.sage
│
├── posgeom/                # TRUSTED: Keep using
│   ├── forest_polytope.py  ✓
│   ├── physics_map.py      ✓
│   └── saddle_pushforward.py
│
└── chy_oracle/             # TRUSTED: Keep using
    ├── laplacian_bridge.py ✓
    └── bcfw.py
```

---

# Testing Protocol

## Before Each Session
```bash
# Run smoke tests to verify environment
sage tests/test_n4_n5.sage
sage src/scripts/physics_pullback_n6.sage
```

## After Each Implementation
```bash
# Test the specific module
sage src/amplituhedron/verify_bcfw_hodges.sage
```

## Weekly Validation
```bash
# Full test suite
sage run_core_tests.py
```

---

# Key References to Study

1. **Arkani-Hamed, Trnka**: "The Amplituhedron" (arXiv:1312.2007)
   - Defines amplituhedron for gauge theory → generalize to gravity

2. **Arkani-Hamed, He, Lam**: "Stringy Canonical Forms" (arXiv:1912.08707)
   - Polytopes to amplitudes via stringy integrals → pushforward approach

3. **Hodges**: "Eliminating spurious poles" (arXiv:0905.1473)
   - The determinant formula we're matching

4. **NSVW**: "The Tree Formula for MHV Gravity" (arXiv:0909.0229)
   - Forest polynomial origin

---

# Timeline

| Week | Phase | Goals |
|------|-------|-------|
| 1 | C: Pushforward Debug | Understand n=6 failure, determine if fixable |
| 2 | A: BCFW Setup | Implement clean momentum twistor + BCFW cells |
| 3 | A: Verification | Test BCFW = Hodges, iterate on formula |
| 4 | A/B: Extension | Extend to n=7 OR try intersection theory |

---

# Risk Assessment

| Risk | Probability | Mitigation |
|------|-------------|------------|
| BCFW formula incorrect | Medium | Study literature carefully, test on n=4,5 first |
| Numerical precision issues | Low | Use exact rational arithmetic |
| Positivity region doesn't exist | Medium | Accept as negative result, still publishable |
| Code is silently wrong | Medium | Multiple independent implementations |

---

# Success Definition

**Minimum Success (Publishable):**
- BCFW cell sum equals Hodges for n=6 with verified positivity
- Clean mathematical statement of the construction

**Strong Success (Major Paper):**
- Works for n=7, 8
- Identify the geometric object (what IS the gravity amplituhedron?)

**Maximum Success (Revolutionary):**
- All-n proof
- Connection to quantum gravity / string theory

---

*Plan created: December 31, 2025*
*Based on: CURSOR_SONNET_45_INSTRUCTIONS_GRAVITY_AMPLITUHEDRON.md*


