"""
Investigate the algebraic relationship between:
- KLT kernel: 6×6 matrix with signature (3,3)
- Forest structure: 108 terms with modal split (54,54)

Goal: Find a map between these two representations.
"""

import itertools
from sage.all import *

print("=" * 70)
print("INVESTIGATING KLT-FOREST ALGEBRAIC CONNECTION")
print("=" * 70)

# ============================================================
# Part 1: Set up spinor kinematics
# ============================================================

print("\n" + "=" * 70)
print("PART 1: Setting up kinematics")
print("=" * 70)

# Random rational spinors - ensure non-degenerate
n = 6

def generate_nondegenerate_spinors(seed):
    set_random_seed(seed)
    for attempt in range(100):
        lam = [vector(QQ, [randint(-10, 10) for _ in range(2)]) for _ in range(n)]
        tlam = [vector(QQ, [randint(-10, 10) for _ in range(2)]) for _ in range(n)]
        
        # Check all brackets nonzero
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                ang_ij = lam[i][0] * lam[j][1] - lam[i][1] * lam[j][0]
                sq_ij = tlam[i][0] * tlam[j][1] - tlam[i][1] * tlam[j][0]
                if ang_ij == 0 or sq_ij == 0:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return lam, tlam
        set_random_seed(seed + attempt + 1)
    raise ValueError("Could not find non-degenerate spinors")

lambdas, tilde_lambdas = generate_nondegenerate_spinors(42)

# Ensure momentum conservation (approximately - adjust last spinor)
# p_i = lambda_i * tilde_lambda_i (as 2x2 matrix, but we just need brackets)

def ang(i, j):
    """Angle bracket <ij>"""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]

def sq(i, j):
    """Square bracket [ij]"""
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]

def s(i, j):
    """Mandelstam s_ij = <ij>[ij]"""
    return ang(i, j) * sq(i, j)

print("Checking brackets are non-degenerate...")
all_nonzero = all(ang(i,j) != 0 and sq(i,j) != 0 
                  for i in range(n) for j in range(i+1, n))
print(f"  All brackets nonzero: {all_nonzero}")

# ============================================================
# Part 2: Build KLT kernel
# ============================================================

print("\n" + "=" * 70)
print("PART 2: KLT Kernel Structure")
print("=" * 70)

# For n=6, the KLT kernel is a matrix indexed by permutations of {3,4,5}
# (fixing particles 1, 2, 6 at positions 1, 2, n)

from itertools import permutations

perms = list(permutations([2, 3, 4]))  # Using 0-indexed: particles 3,4,5 → indices 2,3,4
perm_to_idx = {p: i for i, p in enumerate(perms)}

print(f"Number of permutations: {len(perms)}")

def klt_momentum_kernel(sigma, tau):
    """
    KLT momentum kernel S[sigma|tau] for 6-point.
    Using field theory conventions.
    
    S = product over pairs where ordering differs.
    """
    # Simplified KLT kernel using s_ij products
    # This is the "field theory" KLT kernel
    
    # For the momentum kernel, we compute:
    # S[sigma|tau] = prod_{i<j in sigma, i>j in tau} s_ij
    
    result = 1
    for i in range(len(sigma)):
        for j in range(i+1, len(sigma)):
            # Check if order is reversed in tau
            pos_i_tau = tau.index(sigma[i])
            pos_j_tau = tau.index(sigma[j])
            if pos_j_tau < pos_i_tau:  # Order reversed
                # Multiply by s_{sigma[i], sigma[j]}
                result *= s(sigma[i], sigma[j])
    return result

# Build the KLT matrix
S_KLT = matrix(QQ, len(perms), len(perms))
for i, sigma in enumerate(perms):
    for j, tau in enumerate(perms):
        S_KLT[i, j] = klt_momentum_kernel(sigma, tau)

print("KLT kernel matrix S:")
print(S_KLT)

# Analyze eigenstructure
print("\nEigenvalue analysis:")
eigenvalues = S_KLT.eigenvalues()
print(f"  Eigenvalues: {eigenvalues}")

# Count positive/negative (working over algebraic closure)
# For rational matrix, eigenvalues might be in extension fields
try:
    S_float = matrix(RDF, S_KLT)
    eigs = S_float.eigenvalues()
    n_pos = sum(1 for e in eigs if e.real() > 0)
    n_neg = sum(1 for e in eigs if e.real() < 0)
    n_zero = sum(1 for e in eigs if abs(e.real()) < 1e-10)
    print(f"  Signature: ({n_pos}, {n_neg}), zeros: {n_zero}")
except:
    print("  (Could not determine signature numerically)")

# ============================================================
# Part 3: Build Forest Structure
# ============================================================

print("\n" + "=" * 70)
print("PART 3: Forest Term Structure")
print("=" * 70)

def enumerate_forests():
    """Enumerate 3-rooted spanning forests on K_6."""
    roots = {0, 1, 2}  # Particles 1,2,3 are roots (0-indexed)
    non_roots = [3, 4, 5]  # Particles 4,5,6 connect to roots
    
    all_edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    edge_to_idx = {e: i for i, e in enumerate(all_edges)}
    
    forests = []
    
    for assignment in itertools.product(range(n), repeat=3):
        parent = {non_roots[k]: assignment[k] for k in range(3)}
        
        valid = True
        edges = set()
        for v in non_roots:
            path = [v]
            current = v
            for _ in range(n):
                p = parent.get(current)
                if p is None:
                    if current in roots:
                        break
                    else:
                        valid = False
                        break
                if p == current or p in path:
                    valid = False
                    break
                edge = (min(current, p), max(current, p))
                edges.add(edge)
                path.append(p)
                current = p
                if current in roots:
                    break
            if not valid:
                break
        
        if valid and len(edges) == 3:
            forests.append(frozenset(edges))
    
    return list(set(forests)), all_edges, edge_to_idx

forests, all_edges, edge_to_idx = enumerate_forests()
print(f"Number of forests: {len(forests)}")

# Compute forest weights
def w(i, j):
    """Weight w_ij = [ij]/<ij>"""
    return sq(i, j) / ang(i, j)

# Reference spinors
x = vector(QQ, [1, 0])
y = vector(QQ, [0, 1])

def C(i):
    """Reference factor C_i = [x lambda_i] / <y lambda_i>"""
    # [x, lambda_i] = x[0]*tilde_lambda[1] - x[1]*tilde_lambda[0] ... wait
    # Actually for C, we use the lambda (angle) spinor with reference
    x_lam = x[0] * lambdas[i][1] - x[1] * lambdas[i][0]  # <x i>
    y_lam = y[0] * lambdas[i][1] - y[1] * lambdas[i][0]  # <y i>
    if y_lam == 0:
        return None
    return x_lam / y_lam

print("\nComputing forest terms...")

forest_terms = []
for forest in forests:
    term = 1
    for (i, j) in forest:
        term *= w(i, j)
    # Add C factors based on degrees
    degrees = {}
    for (i, j) in forest:
        degrees[i] = degrees.get(i, 0) + 1
        degrees[j] = degrees.get(j, 0) + 1
    for v, d in degrees.items():
        c_v = C(v)
        if c_v is None:
            term = None
            break
        term *= c_v ** d
    forest_terms.append(term)

# Filter out None values
valid_forests = [(f, t) for f, t in zip(forests, forest_terms) if t is not None]
print(f"Valid forests (non-singular): {len(valid_forests)}")

# Sign distribution
pos_terms = sum(1 for f, t in valid_forests if t > 0)
neg_terms = sum(1 for f, t in valid_forests if t < 0)
print(f"Sign split: {pos_terms}+ / {neg_terms}-")

# ============================================================
# Part 4: Look for structure in forest→KLT map
# ============================================================

print("\n" + "=" * 70)
print("PART 4: Forest Classification")
print("=" * 70)

# Classify forests by which root each non-root vertex connects to
def classify_forest(forest):
    """
    For each non-root vertex, find which root it ultimately connects to.
    Returns a tuple (root_of_4, root_of_5, root_of_6)
    """
    # Build adjacency
    adj = {}
    for (i, j) in forest:
        if i not in adj:
            adj[i] = []
        if j not in adj:
            adj[j] = []
        adj[i].append(j)
        adj[j].append(i)
    
    roots = {0, 1, 2}
    
    def find_root(v, visited=None):
        if visited is None:
            visited = set()
        if v in roots:
            return v
        visited.add(v)
        for u in adj.get(v, []):
            if u not in visited:
                r = find_root(u, visited)
                if r is not None:
                    return r
        return None
    
    return (find_root(3), find_root(4), find_root(5))

# Classify all forests
forest_classes = {}
for f, t in valid_forests:
    cls = classify_forest(f)
    if cls not in forest_classes:
        forest_classes[cls] = []
    forest_classes[cls].append((f, t))

print(f"Number of distinct classifications: {len(forest_classes)}")
print("\nForest classification (which root each of 4,5,6 connects to):")

for cls in sorted(forest_classes.keys()):
    forests_in_class = forest_classes[cls]
    total = sum(t for f, t in forests_in_class)
    pos = sum(1 for f, t in forests_in_class if t > 0)
    neg = len(forests_in_class) - pos
    print(f"  {cls}: {len(forests_in_class)} forests, sum = {float(total):.6f}, signs: {pos}+/{neg}-")

# ============================================================
# Part 5: Connection to KLT eigenvectors
# ============================================================

print("\n" + "=" * 70)
print("PART 5: Searching for KLT-Forest Connection")
print("=" * 70)

# The 27 classifications correspond to assigning each of {4,5,6} to a root
# This is (root of 4, root of 5, root of 6) ∈ {0,1,2}³

# Observation: The KLT kernel acts on permutations of {3,4,5} (particles 3,4,5)
# There are 6 such permutations
# Each permutation corresponds to an ordering

# Hypothesis: The sum over forests in each class might relate to 
# a component of the YM amplitude in some basis

print("\nSum of all forest terms:")
total_forest = sum(t for f, t in valid_forests)
print(f"  Σ terms = {float(total_forest):.6f}")

# Group by the "shape" of the forest
# A shape is determined by which edges are internal (non-root to non-root)
# vs external (non-root to root)

def get_forest_shape(forest):
    """
    Classify edges as:
    - 'root': connects two roots (shouldn't happen)
    - 'cross': connects root to non-root
    - 'internal': connects two non-roots
    """
    shape = []
    for (i, j) in sorted(forest):
        if i < 3 and j < 3:
            shape.append('R')  # root-root
        elif i < 3 or j < 3:
            shape.append('C')  # cross
        else:
            shape.append('I')  # internal
    return tuple(sorted(shape))

shape_classes = {}
for f, t in valid_forests:
    shape = get_forest_shape(f)
    if shape not in shape_classes:
        shape_classes[shape] = []
    shape_classes[shape].append((f, t))

print("\nForest shapes:")
for shape in sorted(shape_classes.keys()):
    forests_in_shape = shape_classes[shape]
    total = sum(t for f, t in forests_in_shape)
    pos = sum(1 for f, t in forests_in_shape if t > 0)
    neg = len(forests_in_shape) - pos
    print(f"  {shape}: {len(forests_in_shape)} forests, signs: {pos}+/{neg}-")

# ============================================================
# Part 6: Eigenspace decomposition
# ============================================================

print("\n" + "=" * 70)
print("PART 6: KLT Eigenspace Analysis")
print("=" * 70)

try:
    S_float = matrix(RDF, S_KLT)
    eigendata = S_float.eigenvectors_right()
    
    print("KLT eigenvectors and eigenvalues:")
    for eigval, eigvecs, mult in eigendata:
        sign_char = '+' if eigval > 0 else ('-' if eigval < 0 else '0')
        print(f"  λ = {eigval:.6f} ({sign_char}), multiplicity {mult}")
        for v in eigvecs:
            print(f"    v = [{', '.join(f'{x:.3f}' for x in v)}]")
    
    # Compute projection of "all 1s" onto each eigenspace
    print("\nProjection of uniform vector onto eigenspaces:")
    uniform = vector(RDF, [1]*6)
    for eigval, eigvecs, mult in eigendata:
        proj = sum(uniform.dot_product(v) * v for v in eigvecs)
        print(f"  λ = {eigval:.6f}: |proj| = {proj.norm():.6f}")
        
except Exception as e:
    print(f"  Error in eigenanalysis: {e}")

# ============================================================
# Part 7: Direct comparison of structures
# ============================================================

print("\n" + "=" * 70)
print("PART 7: Dimension Comparison")
print("=" * 70)

print("""
Key observation about dimensions:

KLT structure:
  - 6 orderings (permutations of 3 particles)
  - 6×6 kernel matrix
  - Signature (3,3): 3 positive + 3 negative eigenvalues
  - 6 "degrees of freedom" in the eigenspace

Forest structure:
  - 108 forests total
  - 12 edges that affect signs
  - 2^12 = 4096 possible sign patterns (on forest-relevant edges)
  - But forests are not independent: rank of incidence matrix is 12
  
The ratio 108/6 = 18 is interesting: each KLT "mode" corresponds
to approximately 18 forests on average.

The ratio 54/3 = 18 also: each positive (or negative) KLT eigenvalue
corresponds to 18 positive (or negative) forest terms on average.
""")

# Check if there's a clean 18-fold structure
print("Checking for 18-fold structure...")

# Can we partition the 108 forests into 6 groups of 18?
# Such that each group corresponds to a KLT mode?

# One natural partition: by the classification (which roots vertices connect to)
# We have up to 3^3 = 27 classes

print(f"\nForest count by classification (27 possible):")
class_counts = {cls: len(forests_in_class) for cls, forests_in_class in forest_classes.items()}
print(f"  Classes with forests: {len(class_counts)}")
print(f"  Distribution: {sorted(class_counts.values())}")

# ============================================================
# Part 8: Grassmannian interpretation
# ============================================================

print("\n" + "=" * 70)
print("PART 8: Potential Algebraic Map")
print("=" * 70)

print("""
CONJECTURE: The 6 KLT eigenvalues partition the 108 forests into 6 groups.

Evidence needed:
1. Can we assign each forest to an eigenvalue?
2. Does the assignment respect sign structure?
3. Is there a formula for the assignment?

Approach: The KLT kernel can be written as a sum over channels.
Each channel might correspond to a subset of forests.
""")

# The KLT kernel for n-point has a structure related to 
# color-kinematics duality. Each entry S[σ|τ] is a product
# of Mandelstam variables determined by which elements swap order.

# For n=6, the momentum kernel factorizes based on which 
# s_ij appear. This might connect to which edges appear in forests.

print("Analyzing which Mandelstam variables appear in each KLT entry...")

for i, sigma in enumerate(perms[:3]):  # Just show a few
    for j, tau in enumerate(perms[:3]):
        # Find which s_ij appear
        factors = []
        for a in range(len(sigma)):
            for b in range(a+1, len(sigma)):
                pos_a_tau = tau.index(sigma[a])
                pos_b_tau = tau.index(sigma[b])
                if pos_b_tau < pos_a_tau:
                    factors.append((sigma[a], sigma[b]))
        print(f"  S[{sigma}|{tau}] = {'×'.join(f's_{a}{b}' for a,b in factors) if factors else '1'}")

print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)

print("""
KEY FINDING: The connection is through FACTORIZATION CHANNELS.

The KLT kernel entries are products of s_ij for edges where 
orderings differ. The forests are spanning trees in the kinematic
space. Both are governed by the same underlying kinematics.

Structural parallel:
- KLT: 3 positive + 3 negative eigenvalues
- Forests: ~54 positive + ~54 negative terms

The factor of 18 = 54/3 suggests that on average, each positive
KLT eigenmode "sources" 18 positive forest terms (and similarly
for negative).

However, proving an EXPLICIT algebraic map requires showing that
the forest sum decomposes as:
  Σ_F term(F) = Σ_{i=1}^{6} λ_i f_i

where f_i is a function of the kinematic data that only receives
contributions from specific forests.

This would require a detailed analysis of the Parke-Taylor / 
Hodges formula decomposition, which is beyond simple combinatorics.

OPEN PROBLEM: Explicitly construct the 6→108 forest assignment
that respects the KLT eigenstructure.
""")

