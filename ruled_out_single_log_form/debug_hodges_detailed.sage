load('proper_proof_amplituhedron_hodges.sage')

Z = sample_positive_Z_moment_curve(n=6, seed=0)
twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)

print("="*70)
print("DETAILED HODGES DEBUG")
print("="*70)

# Check all angle brackets used in Hodges
print("\nAngle brackets for off-diagonal Phi:")
for i in range(6):
    for j in range(6):
        if i != j:
            ang = twistor.get_angle(i, j)
            if ang == 0:
                print(f"  ZERO: <{i},{j}> = {ang}")

print("\nReference leg brackets (x=0, y=5):")
for i in range(6):
    ang_i0 = twistor.get_angle(i, 0)
    ang_i5 = twistor.get_angle(i, 5)
    print(f"  <{i},0> = {ang_i0}, <{i},5> = {ang_i5}")
    if ang_i0 == 0 or ang_i5 == 0:
        print(f"    WARNING: Zero bracket for particle {i}")

print("\nNormalization brackets:")
ang_01 = twistor.get_angle(0, 1)
ang_12 = twistor.get_angle(1, 2)
ang_20 = twistor.get_angle(2, 0)
print(f"  <0,1> = {ang_01}")
print(f"  <1,2> = {ang_12}")
print(f"  <2,0> = {ang_20}")
if ang_01 == 0 or ang_12 == 0 or ang_20 == 0:
    print("    WARNING: Zero normalization bracket")

print("\nConsecutive brackets:")
for i in range(6):
    j = (i + 1) % 6
    ang = twistor.get_angle(i, j)
    print(f"  <{i},{j}> = {ang}")
    if ang == 0:
        print(f"    WARNING: Zero consecutive bracket")

# Try computing Hodges step by step
print("\n" + "="*70)
print("STEP-BY-STEP HODGES COMPUTATION")
print("="*70)

# Build Phi matrix
n = 6
Phi = matrix(QQ, n, n)
x, y = 0, 5

# Off-diagonal
print("\nBuilding off-diagonal Phi...")
for i in range(n):
    for j in range(n):
        if i != j:
            ij_ang = twistor.get_angle(i, j)
            if ij_ang == 0:
                print(f"  ERROR: <{i},{j}> = 0")
            ij_sq = twistor.get_square(i, j)
            if ij_sq is None:
                print(f"  ERROR: [i,j] undefined for ({i},{j})")
            else:
                Phi[i, j] = ij_sq / ij_ang

print("Off-diagonal complete")

# Diagonal
print("\nBuilding diagonal Phi...")
for i in range(n):
    ix_ang = twistor.get_angle(i, x)
    iy_ang = twistor.get_angle(i, y)
    if ix_ang == 0 or iy_ang == 0:
        print(f"  ERROR: Reference brackets zero for particle {i}")
    else:
        diag_sum = QQ(0)
        for j in range(n):
            if j == i:
                continue
            jx_ang = twistor.get_angle(j, x)
            jy_ang = twistor.get_angle(j, y)
            if jx_ang == 0 or jy_ang == 0:
                continue
            contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
            diag_sum -= contrib
        Phi[i, i] = diag_sum

print("Diagonal complete")
print(f"Phi matrix built: {Phi.nrows()}x{Phi.ncols()}")

# Reduced determinant
print("\nComputing reduced determinant...")
rows_to_keep = [3, 4, 5]
cols_to_keep = [3, 4, 5]
Phi_red = Phi[rows_to_keep, cols_to_keep]
det_Phi_red = Phi_red.det()
print(f"det(Phi_red) = {det_Phi_red}")

# Normalization
norm_factor = (ang_01 * ang_12 * ang_20) ** 2
print(f"norm_factor = {norm_factor}")

if norm_factor == 0:
    print("  ERROR: Normalization factor is zero")
else:
    det_prime_Phi = det_Phi_red / norm_factor
    print(f"det'(Phi) = {det_prime_Phi}")

