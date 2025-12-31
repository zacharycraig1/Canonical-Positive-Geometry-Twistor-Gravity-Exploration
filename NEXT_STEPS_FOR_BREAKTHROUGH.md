# Next Steps for Positive Geometry Gravity Breakthrough

## Current State
We have fixed the foundational linear algebra objects:
- Hodges Matrix $\Phi$ has correct rank 3.
- KLT Kernel $S$ has correct rank 6.
- Both amplitudes scale correctly ($Z^{-8}$).
- Ratio is $\mathcal{O}(1)$.

## Roadmap

### Step 1: Pin Down the Exact Constant
The ratio varies (1.3 - 1.6). We must make it **exactly 1** (or a constant).
- **Action**: Check if the variation correlates with a cross-ratio $u = \frac{s_{12}s_{45}}{s_{14}s_{25}}$.
- **Action**: Verify the Parke-Taylor normalization factors in KLT sum.
- **Goal**: `assert(abs(ratio - 1.0) < 1e-10)`

### Step 2: The Intersection Matrix (The "Metric")
In intersection theory, the KLT kernel $S$ acts as the metric between "Parke-Taylor forms".
$$ \langle PT(\alpha) | PT(\beta) \rangle = S[\alpha|\beta]^{-1} $$
- **Action**: Compute the matrix of "Intersection Numbers" for logarithmic forms on the moduli space $\mathcal{M}_{0,6}$.
- **Breakthrough**: Show that $S_{KLT}$ is exactly the inverse of the Intersection Matrix of the Associahedron canonical forms.

### Step 3: The "Squared" Geometry
Gravity is "YM squared".
- **Concept**: The geometry for gravity is not a simple polytope in kinematic space.
- **Hypothesis**: It is a "Weighted Positive Geometry".
- **Action**: Construct the "Weighted Canonical Form":
  $$ \Omega_{Grav} = \sum_{\alpha, \beta} S_{\alpha \beta} \Omega_\alpha \Omega_\beta $$
  where $\Omega_\alpha$ is the canonical form of the permuted Associahedron.

### Step 4: Geometric Realization (The "Amplituhedron" for Gravity)
Can we find a *single* geometric object $\mathcal{G}$?
- Recent literature suggests **Cayley Polytopes** or **Hessian Geometry**.
- **Action**: Map the KLT kernel to the "Hessian of the volume" of the Associahedron.

## Immediate Command
Run the `check_normalization` search to fix the O(1) drift.
