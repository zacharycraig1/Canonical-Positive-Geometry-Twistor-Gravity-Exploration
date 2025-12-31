# Analysis of Hodges Numerator and Poles

## 1. Pole Order at $\langle i, i+1 \rangle \to 0$
- **Result:** Estimated Pole Order is **-6.00**.
- **Interpretation:** The amplitude scales as $M \sim (\langle i, i+1 \rangle)^{-6}$ near the pole?
  - Wait, standard gravity MHV usually scales as $1/t^3$ or $1/s^3$?
  - For MHV, $M_{YM} \sim 1/\langle 1 2 \rangle$. $M_{Gr} \sim M_{YM}^2 \times S$.
  - This suggests effective pole order is high.
  - The denominator used in my previous ansatz was $(\prod \langle i, i+1 \rangle)^2$. Each factor contributes order 2.
  - If pole order is 6, it means $M \sim 1/\langle 01 \rangle^6$?
  - But $\prod^2$ would give $1/\langle 01 \rangle^2$.
  - So the Numerator $N$ must scale as $\langle 01 \rangle^{-4}$? This contradicts polynomiality.
  - **Correction:** The "Estimated Pole Order" calculation was naive (slope of log).
    - If $M \propto 1/t^k$, slope is $-k$.
    - Result -6.00 means $M \propto 1/t^6$.
    - The denominator factor $(\prod \langle i i+1 \rangle)^2$ contributes $1/t^2$ (since only one factor vanishes).
    - So $N = M \times D$ scales as $1/t^6 \times t^2 = 1/t^4$.
    - This suggests **N is singular** at collinear limits? This contradicts polynomiality.
    - **Hypothesis:** The "Soft Limit" involves more than one bracket vanishing, or the "Pole Order" check was contaminated by soft limit behavior ($\lambda \to 0$ scales all brackets).
    - Here I perturbed $Z[1]$ towards $Z[0]$. This implies $\langle 0 1 \rangle \to 0$.
    - But it also implies $\langle 1 2 \rangle \to \langle 0 2 \rangle$ (finite).
    - So only one bracket in the cyclic product vanishes.
    - If $M \sim 1/\langle 0 1 \rangle^6$, then gravity is very singular.
    - Standard result: Gravity amplitudes have $1/p^2$ poles. In spinor variables, $s_{01} = \langle 01 \rangle [01]$.
    - If $[01]$ is finite, then pole is $1/\langle 01 \rangle$.
    - But BCFW shifts usually show $1/z$ poles.
    - **Re-evaluate:** Maybe the denominator is $(\prod \langle i, i+1 \rangle)^6$?? Unlikely.

## 2. Spurious Poles $\langle 0, 2 \rangle \to 0$
- **Result:** Estimated Spurious Pole Order is **-0.01**.
- **Interpretation:** The amplitude is **finite** (order 0) as $\langle 0, 2 \rangle \to 0$.
- **Conclusion:** There are **NO spurious poles**. This is consistent with physical expectations.
- **Koefler Condition:** "Numerator vanishes whenever any two spinor brackets go to zero".
  - My check showed $N \approx -2.5e5$ (non-zero) at $\langle 0,2 \rangle \to 0$.
  - Wait, if M is finite and $D = (\prod \langle i i+1 \rangle)^2$ is finite (since $\langle 0,2 \rangle$ is not in the product), then $N = M \times D$ must be finite (non-zero).
  - So $N$ does **NOT** vanish at $\langle 0,2 \rangle \to 0$.
  - This contradicts the prompt's interpretation of Koefler et al.: "enforce that N vanishes whenever any two spinor brackets go to zero".
  - **Correction:** Perhaps the Koefler condition refers to a *different* denominator?
  - Koefler et al. (arXiv:1905.00849? No, maybe 20xx) "Gravity MHV ... numerator is a nontrivial polynomial with specified zeros".
  - If the denominator was $\prod_{all} \langle i j \rangle$, then N would need to vanish to cancel spurious poles.
  - But we found M is finite at spurious locations.
  - So if we use the *minimal* denominator $D_{cyclic}$, N need not vanish.
  - BUT, does N vanish at physical poles to reduce the order?
  - We found $M \sim 1/t^6$ at physical pole. $D_{cyclic} \sim t^2$.
  - This implies $N \sim 1/t^4$. This is a huge problem for the "Polynomial Numerator" claim.

## 3. Resolving the Pole Order Discrepancy
- The "Polynomial N" result from the previous step (Degree 34) was very strong.
- How can N be a polynomial if it blows up at $\langle 01 \rangle \to 0$?
- **Possibility:** The denominator $D_{cyclic}$ is NOT just $(\prod \langle i i+1 \rangle)^2$.
- The scaling weight of $D_{cyclic}$ is $6 \times 2 = 12 \times 2 = 24$ (in lambda?). No, $\langle i j \rangle$ is weight 1 in $Z_i$. So weight 2 per bracket. Total weight 24.
- $M_{MHV}$ has weight -2 (Dimension Z^-2).
- So N must have weight 22.
- The degree 34 polynomial in $t$ (where $Z(t)$ is linear) matches this?
  - $N(Z)$ is homogeneous of degree $k$. $N(Z_A + t Z_B)$ is degree $k$.
  - My previous script found degree 34.
  - This implies weight 34?
  - Discrepancy: Weight 22 vs Degree 34.
  - **Explanation:** $Z_i$ has 4 components. Degree in $t$ depends on how $Z(t)$ is constructed.
  - I used $Z_i(t) = Z0_i + t(Z1_i - Z0_i)$. This is linear in $t$.
  - So degree in $t$ = Homogeneous degree in $Z$.
  - So $N$ has weight 34.
  - Denominator $(\prod \langle i i+1 \rangle)^2$ has weight $12 \times 2 = 24$?
  - $\langle i j \rangle$ is bilinear in Z? No, quadratic in Z?
  - $\langle i j \rangle = Z_i \cdot I \cdot Z_j$ ? No.
  - $\langle i j \rangle$ is determinant of 2x2 matrix of top components. It is bilinear in $Z$.
  - So weight 2.
  - $(\prod_{i=1}^6 \langle i i+1 \rangle)^2$ has weight $6 \times 1 \times 2 = 12$?
  - No, product of 6 brackets. Each bracket is deg 1 in $Z_i$.
  - Total degree in all Z is $6 \times 2 = 12$.
  - Squared: 24.
  - So Weight(D) = 24.
  - If Weight(N) = 34.
  - Then Weight(M) = 34 - 24 = 10.
  - But Gravity Amplitude must have weight -2 (in twistors, usually scaling Z scales M?).
  - In spinor helicity, $M(\lambda, \tilde\lambda)$ scales as $\lambda^{-2}$?
  - Twist: $Z \to c Z$. M scales as $c^k$.
  - My scaling check found $k = -2$ (or something close? "Scaling Dimension (log2): -2.00").
  - So M has weight -2.
  - This requires Weight(N) - Weight(D) = -2.
  - Weight(N) = 24 - 2 = 22.
  - **Contradiction:** Polynomial check found Degree 34.
  - **Wait:** In the polynomial check, I interpolated $N(t)$.
  - Did I use the *same* denominator function?
  - `def get_denom(tw): ... return val**2`.
  - Yes.
  - Why did `poly.degree()` return 34?
  - Let's re-read the log. `Found Polynomial of degree 34`.
  - If Weight is 22, why degree 34?
  - Maybe the "Scaling Dimension" check measures something else?
  - Scaling Z -> 2Z.
  - `H2/H1 = 0.25 = 2^-2`. So weight -2. Correct.
  - `D2/D1`: D scales as $Z^{24}$. So $2^{24}$.
  - `N2/N1 = (H2*D2) / (H1*D1) = 2^{-2} * 2^{24} = 2^{22}$.
  - So N SHOULD have degree 22.
  - Why did Lagrange give 34?
  - **Hypothesis:** $N$ is NOT a polynomial in Z, but `check_numerator_structure.sage` found a polynomial fit because I sampled on a line?
  - But I checked extrapolation. It was perfect.
  - A rational function restricts to a rational function on a line. If it matched a polynomial on a line, it IS a polynomial (or denominator is constant).
  - Maybe I am miscounting the degree of D?
  - `val *= tw.get_angle(...)`. Angle is quadratic in Z? No, bilinear. Deg 2 in Z (total).
  - 6 angles => Deg 12.
  - Squared => Deg 24.
  - So Degree 22 is correct.
  - **Why 34?**
  - Maybe `hodges_6pt_mhv` includes a factor I missed?
  - It returns `det_prime * <0 1>^8`.
  - `<0 1>^8` has degree 8 (in angle) = 16 in Z.
  - `det_prime` has weight?
  - `Phi` elements are $[ij]/<ij>$. Weight 0?
  - No. $[ij]$ is weight 2 in dual space? In twistors, $[ij]$ involves derivative or integral?
  - $[ij]$ is defined via dual twistors or derivative.
  - In `src/hodges.sage`, `get_square` uses `get_tilde_lambda`.
  - `get_tilde_lambda` uses `mu` (deg 1) and `angle` (deg 2).
  - Formula: `(mu <...> + ...) / (<><>)`.
  - Num: Deg 1 + Deg 2 = 3. Denom: Deg 4.
  - So `tilde_lambda` is Deg -1.
  - `get_square`: `lam_t_i * lam_t_j`. Deg -2.
  - `get_angle`: Deg 2.
  - `Phi_ij = [ij]/<ij>`. Deg -2 - 2 = -4.
  - `Phi` matrix is 6x6. Rank 3. Reduced det is 3x3.
  - `det(Phi_red)`: 3 * (-4) = -12.
  - `norm_factor`: `(<><><>)^2`. Deg (2+2+2)*2 = 12.
  - `det_prime = det / norm`. Deg -12 - 12 = -24.
  - Total M = `det_prime * <0 1>^8`.
  - Deg = -24 + (8*2) = -24 + 16 = -8.
  - **Wait:** Scaling check said -2.
  - Discrepancy: -8 vs -2.
  - Let's check `get_tilde_lambda` scaling.
  - Z -> cZ. mu -> c mu. angle -> c^2 angle.
  - Num -> c * c^2 = c^3. Denom -> c^4.
  - tilde_lambda -> c^-1. Correct.
  - square -> c^-2. Correct.
  - angle -> c^2. Correct.
  - Phi -> c^-2 / c^2 = c^-4. Correct.
  - det_red (3x3) -> (c^-4)^3 = c^-12. Correct.
  - norm -> (c^6)^2 = c^12. Correct.
  - det_prime -> c^-12 / c^12 = c^-24. Correct.
  - factor <01>^8 -> (c^2)^8 = c^16.
  - Total -> c^-24 * c^16 = c^-8.
  - **Why did scaling check say -2?**
  - "Scaling Dimension (log2): -2.00"
  - Let's look at `src/hodges.sage`:
  - `Phi[i, i] = - sum Phi[ij] ...`
  - Elements of Phi scale as c^-4.
  - So det scales as c^-12.
  - Is `det_Phi_red` really 3x3?
  - n=6. Deleted 3. 3x3. Yes.
  - Maybe I am wrong about `norm_factor`?
  - `norm_factor = (ang*ang*ang)**2`.
  - Maybe `det_prime` scaling is different.
  - **Crucial:** The "Scaling Dimension" output -2.00 is experimental fact.
  - My derivation -8 is theoretical.
  - Where is the factor of c^6?
  - Ah, `<0 1>^8`. Maybe it is `<0 1>^4`?
  - Standard MHV is `<ij>^4`.
  - If factor is `<0 1>^4`, then c^8.
  - -24 + 8 = -16. Still not -2.
  - Maybe Phi elements scale as c^0?
  - $[ij]$ scales as c^2? No.
  - **Twistor Scaling:**
  - $Z \to cZ$.
  - This is a scaling of the projective coordinates.
  - Amplitude should be Weight 0 in projective weights?
  - No, Amplitude has dimension mass^(4-n)?
  - In twistors, usually we fix scaling.
  - But if we check homogeneity, we find the weight.
  - If exp says -2, then N has weight 22 (given D has 24).
  - Why did polynomial fit give 34?
  - 34 - 24 = 10.
  - If M has weight 10.
  - 10 is not -2.
  - **Discrepancy:** The polynomial check `poly.degree()` = 34.
  - But `check_scaling_poly.sage` found `log2(ratio) = -2`.
  - These two facts from the SAME script run are contradictory if N=M*D.
  - `N2/N1 = 2^34` implies `log2(M2*D2 / M1*D1) = 34`.
  - `log2(M2/M1) + log2(D2/D1) = 34`.
  - `-2 + 24 = 22`.
  - So N *should* be degree 22.
  - But `poly.degree()` says 34.
  - **Re-read log:** `Found Polynomial of degree 34`.
  - `Extrapolation to t=40 ... SUCCESS`.
  - This means the degree IS 34.
  - So `log2(N2/N1)` MUST be 34.
  - So `log2(M2/M1)` MUST be 34 - 24 = 10.
  - But the script printed `Scaling Dimension (log2): -2.00`.
  - **Resolution:** I must have used *different* functions or setups for the two checks in `check_scaling_poly.sage`?
  - No, same `hodges_6pt_mhv`.
  - **Wait:** In scaling check, I used `Z2 = 2*Z`.
  - In Poly check, I used `Z(t) = Z0 + t(Z1-Z0)`.
  - Degree in $t$ corresponds to weight in Z ONLY if $N(Z)$ is homogeneous.
  - Is $N(Z)$ homogeneous?
  - Yes, twistors are projective.
  - **Possibility:** `hodges_6pt_mhv` is NOT homogeneous?
  - It depends on `get_square`.
  - `get_square` uses `get_tilde_lambda`.
  - `get_tilde_lambda` depends on `mu` and `angle`.
  - If I scale Z -> cZ, then mu -> c mu, angle -> c^2 angle.
  - Formula: `(mu <...> + ...) / (<><>)`.
  - c^3 / c^4 = c^-1.
  - Correct.
  - **MYSTERY:** How can Scaling be -2 and Degree be 34?
  - 34 - 24 = 10.
  - Is it possible `get_denom` is NOT weight 24?
  - `(Prod <i,i+1>)^2`. 6 brackets. Each weight 2. Total 12. Square 24.
  - This seems robust.
  - Is it possible `M` scaling is actually 10?
  - Let's check the log output again.
  - `H1: -2.7080e+01`
  - `H2: -6.7700e+00`
  - `Ratio: 2.5000e-01`
  - 0.25 = 1/4 = 2^-2.
  - So scaling is definitely -2.
  - So N weight must be 22.
  - Why is Polynomial Degree 34?
  - Maybe $N(Z)$ is *not* homogeneous? But M is homogeneous (-2). D is homogeneous (24). N=MD must be homogeneous (22).
  - **Wait:** $Z(t) = Z0 + t(Z1-Z0)$ is NOT a scaling.
  - It's a line.
  - If $N(Z)$ is homogeneous of degree $k$, then $N(Z0 + t(Z1-Z0))$ is a polynomial in $t$ of degree $k$.
  - Unless $Z1-Z0$ is special? No, random.
  - **ERROR source?**
  - Maybe `hodges_6pt_mhv` uses `reference_spinors`?
  - `sample_reference_spinors`.
  - In `check_scaling_poly`, I called `hodges_6pt_mhv(tw)`.
  - Inside, it calls `sample_reference_spinors(tw)`.
  - `sample_reference_spinors` uses `random`.
  - Does it depend on Z?
  - It generates random spinors `lambda_x, lambda_y` and checks they are valid.
  - **CRITICAL:** The reference spinors are CONSTANT (or random independent of Z scaling)?
  - If I scale Z -> 2Z, do I scale reference spinors?
  - No, they are auxiliary.
  - But Hodges formula depends on them?
  - **Independence:** Hodges formula is independent of reference spinors.
  - So scaling Z should scale M.
  - But does the *polynomial representation* depend on reference spinors?
  - No, if M is unique.
  - **Wait!** The log said `Found Polynomial of degree 34`.
  - If N is weight 22, degree should be 22.
  - Why 34?
  - 34 - 22 = 12.
  - Where does 12 come from?
  - `norm_factor` in Hodges? It has weight 12.
  - Maybe `norm_factor` is somehow cancelling something?
  - **Hypothesis:** `analyze_hodges_structure.sage` pole analysis suggests $M \sim 1/t^6$.
  - $D \sim t^2$.
  - So $N \sim 1/t^4$.
  - This contradicts polynomiality.
  - UNLESS `get_denom` is wrong.
  - Maybe I should use a denominator that cancels the poles?
  - If pole order is 6, we need $D \sim t^6$.
  - So $D = (\prod \langle i i+1 \rangle)^6$.
  - Weight of D = 72.
  - Weight of N = 72 - 2 = 70.
  - Degree 34 is not 70.

  **Conclusion:** I am confused about the degree and pole order.
  But the **fact** remains: $N = M \times (\prod \langle i i+1 \rangle)^2$ fits a polynomial of degree 34.
  The Extrapolation was perfect.
  So $N$ **IS** a polynomial of degree 34.
  And $M$ **IS** scaling as weight -2.
  This implies $D$ must scale as weight 36.
  $34 - 36 = -2$.
  But $(\prod \langle i i+1 \rangle)^2$ scales as 24.
  **Difference:** 12.
  **Re-check scaling:**
  `H2/H1 = 0.25`.
  `D2/D1`: `get_denom` uses `val**2`.
  `val` is product of 6 angles.
  Scale Z -> 2Z. Angle -> 4 * Angle.
  `val` -> 4^6 * val = 4096 * val.
  `val**2` -> (4096)^2 * val^2.
  4096 = 2^12.
  Squared = 2^24.
  So D scales as 2^24.
  M scales as 2^-2.
  N = M*D scales as 2^22.
  **This is mathematically robust.**
  
  **So why Degree 34?**
  Is it possible `hodges_6pt_mhv` output depends on `seed` which changes with `t`?
  No, `seed` was fixed for Z generation.
  Inside `hodges`, `sample_reference_spinors` uses `random.seed(42)`. Constant.
  
  **Possibility:** My polynomial fit is fitting **noise**?
  No, extrapolation was 1.000000.
  
  **Possibility:** $N(Z)$ is a polynomial of degree 34, but implies M has weight +10.
  But scaling check says weight -2.
  
  **THERE MUST BE A BUG IN MY SCALING OR DEGREE LOGIC.**
  
  Let's look at `check_scaling_poly.sage` again.
  `Z2 = [z * QQ(2) for z in Z]`
  `tw2 = MomentumTwistor(n=6, Z=Z2)`
  `H2 = ...`
  This is correct.
  
  `poly.degree()`.
  The points $(t, N)$ lie on a polynomial.
  $Z(t) = Z0 + t \Delta Z$.
  $N(Z(t))$ is a polynomial in $t$.
  If $N(Z)$ is homogeneous of degree $k$, $N(Z(t))$ is degree $k$.
  
  **Is it possible that <0 1>^8 factor is doing something?**
  Scale Z -> 2Z. <0 1> -> 4 <0 1>.
  ^8 -> 4^8 = 2^16.
  Scaling -2.
  
  Maybe `det_prime` scales as 2^(-18)?
  -18 + 16 = -2.
  
  This consistency check is driving me crazy but it is essential.
  Wait, I see `norm_factor = (ang*ang*ang)**2`.
  Weight 12.
  Maybe I calculated weight of `det_Phi_red` wrong.
  Phi elements scale as $c^{-4}$.
  3x3 det scales as $c^{-12}$.
  $c^{-12} / c^{12} = c^{-24}$.
  Correct.
  
  **Maybe I am counting the degree of $t$ wrong?**
  $Z(t)$ is linear in $t$.
  $\langle i j \rangle (t)$ is bilinear in $Z(t)$. So quadratic in $t$.
  Phi elements $[ij]/<ij>$.
  Num $[ij]$: deg -2 in Z?
  No, let's trace $t$ degree.
  Z -> $t^1$.
  angle -> $t^2$.
  mu -> $t^1$.
  tilde_lambda -> $t^1 * t^2 / t^2 = t^1$?
  Let's recheck `get_tilde_lambda`.
  Num: `mu * angle`. $t^1 * t^2 = t^3$.
  Denom: `angle * angle`. $t^4$.
  Result: $t^{-1}$.
  Square $[ij]$: $t^{-1} * t^{-1} = t^{-2}$.
  Phi element: $[ij]/<ij> = t^{-2} / t^2 = t^{-4}$.
  Det (3x3): $(t^{-4})^3 = t^{-12}$.
  Norm factor: $(t^2 * t^2 * t^2)^2 = t^{12}$.
  Reduced Det: $t^{-12} / t^{12} = t^{-24}$.
  Factor $<01>^8$: $(t^2)^8 = t^{16}$.
  Total M: $t^{-8}$.
  
  **So M scales as $t^{-8}$.**
  **Scaling check said $t^{-2}$.**
  
  WHY?
  Let's look at scaling check again.
  `ratio_H = H2 / H1`.
  `log2(ratio_H)`.
  If $H \sim t^{-8}$, ratio should be $2^{-8} = 1/256 \approx 0.0039$.
  Log should be -8.
  Output was -2.
  
  **What is wrong?**
  Maybe `Phi[i,i]` diagonal calculation?
  `diag_sum -= Phi[ij] * (ang * ang) / (ang * ang)`.
  The fraction cancels weights.
  So diagonal has same weight as off-diagonal. $t^{-4}$.
  
  Maybe `get_tilde_lambda` scaling?
  `num = mu_im1 * ang + ...`
  If Z -> cZ.
  mu -> c mu. ang -> c^2 ang.
  num -> c^3.
  denom -> c^4.
  result -> c^-1.
  This seems solid.
  
  **Is it possible MomentumTwistor normalization changes?**
  No.
  
  **Maybe I am calculating <0 1> scaling wrong?**
  $\langle i j \rangle = Z_i \cdot I \cdot Z_j$.
  $Z_i = (z0, z1, z2, z3)$.
  $<i j> = z0_i z1_j - z1_i z0_j$.
  If $Z \to c Z$, then $<i j> \to c^2 <i j>$.
  Correct.
  
  **Maybe the factor is NOT $<0 1>^8$?**
  In `hodges.sage`: `helicity_factor = twistor.get_angle(0, 1) ** 8`.
  
  **Wait, if scaling is -2, then N weight is 22. But poly degree is 34.**
  **34 = 22 + 12.**
  **12 is the weight of `norm_factor`.**
  Maybe I am NOT dividing by `norm_factor` in my polynomial check?
  No, I call `hodges_6pt_mhv` which does the division.
  
  **Let's look at the Pole Order analysis.**
  Order 6.
  This means $M \sim 1/t^6$ near pole?
  The pole is $\langle i i+1 \rangle \to 0$.
  In "t-space" (Z perturbation), $\langle i i+1 \rangle$ goes as $t$ (linear)?
  No, if $Z1 = Z0 + t(Z1-Z0)$, then $\langle 0 1 \rangle \sim t$?
  Let's check.
  $Z0$ is fixed. $Z1(t) \to Z0$ as $t \to 0$ (if $Z1(0)=Z0$).
  My perturbation: $Z_mod[1] = Z[0] + t*(Z[1]-Z[0])$.
  So $Z_mod[1] \to Z[0]$.
  $\langle 0, 1_{mod} \rangle = \det(Z0, Z0+t\Delta) = \det(Z0, t\Delta) = t \det(Z0, \Delta) = t \langle 0, 1_{orig} \rangle$.
  So bracket is linear in t.
  So Pole Order 6 means $M \sim 1/t^6$.
  $D \sim t^2$ (from $(\prod)^2$, only one vanishes).
  $N = M D \sim t^{-4}$.
  Singular.
  
  **BUT Extrapolation was Perfect.**
  A singular function cannot be extrapolated by a polynomial perfectly.
  This implies $N$ is **NOT** singular.
  This implies my Pole Order estimate is wrong or $N$ cancels the singularity.
  
  **Hypothesis:** The "Pole Order 6" is an artifact of **Collinear Limits**.
  When $Z1 \to Z0$, $\langle 0 1 \rangle \to 0$.
  But also $\langle 1 2 \rangle \to \langle 0 2 \rangle$.
  Maybe other things vanish?
  In collinear limit $1 \to 0$, $s_{01} \to 0$.
  Gravity amplitudes have $1/s^3$ singularities?
  Usually $M \sim 1/\langle 01 \rangle [01]$.
  If $[01] \to 0$ too?
  $[01]$ depends on derivatives.
  If $Z1 \to Z0$, then curve is tangent?
  
  **Okay, I need to trust the Polynomial Fit.**
  It worked on a random line.
  Random lines do not hit singularities.
  So generically $N$ is a polynomial.
  The degree 34 is the generic degree.
  
  **Task: Construct N explicitly.**
  I cannot fit 34 degree polynomial.
  But I can use the "Uniqueness Constraints" again.
  The user said: "set up the ansatz that the 6-pt amplitude equals some numerator N over the known denominators".
  "enforce that N vanishes..."
  
  Maybe I should fit the ratio $M_{KLT}/M_{Hodges}$ again, but now correcting for the "missing polynomial".
  I found `Ratio ~ 1.02` (varies).
  If $M_{KLT} = N_{KLT}/D$.
  $M_{Hodges} = N_{Hodges}/D$.
  $Ratio = N_{KLT} / N_{Hodges}$.
  If N_Hodges is degree 34, maybe N_KLT is also degree 34?
  KLT terms are products of PT factors and s_ij.
  PT factor: $1/\langle\rangle^6$.
  s_ij: $\langle\rangle[]$.
  KLT kernel: $s^3$.
  Total: $1/\langle\rangle^6 * 1/\langle\rangle^6 * s^3 \sim \langle\rangle^{-12} * \langle\rangle^3 []^3$.
  This looks like rational function, not polynomial numerator.
  
  **Wait.**
  The "Polynomial Numerator" is a property of the **reduced** amplitude over a common denominator.
  $M_{KLT}$ can be brought to common denominator $D_{cyclic}$.
  Then its numerator $N_{KLT}$ would be a polynomial.
  Difference $N_{Diff} = N_{Hodges} - N_{KLT}$.
  This $N_{Diff}$ is what we need to find.
  
  **Strategy:**
  1. Calculate $N_{KLT} = M_{KLT} \times D_{cyclic}$.
  2. Calculate $N_{Hodges} = M_{Hodges} \times D_{cyclic}$.
  3. Calculate $\Delta = N_{Hodges} - N_{KLT}$.
  4. Check if $\Delta$ is small or has simple structure.
     - Is $\Delta$ a polynomial?
     - What is its degree?
     - Does it factorize?

  I will write `analyze_numerator_diff.sage`.
  It will sample $\Delta$ on a line and check its degree.
  If $\Delta$ is simpler than $N_{Hodges}$ (lower degree?), that's a hint.
  
  Also, I need to check the KLT Kernel implementation.
  The user asked to "Refine the KLT Double-Copy Implementation".
  Maybe I should try the "Associahedral Grid" kernel or similar.
  But sticking to the standard kernel first.
  Is there a "sign" error?
  The prompt says "check sign conventions".
  I checked permutation symmetry and it failed.
  Maybe I should try to **restore symmetry** by summing over permutations?
  My `check_full_symmetrization` failed (ratio was huge).
  This implies the "standard" basis terms are huge and cancel out carefully?
  KLT usually involves minimal basis $(n-3)!$.
  For n=6, $3! = 6$ terms.
  My code summed 36 terms (6 alpha, 6 beta).
  This is correct for the symmetric sum form:
  $M = \sum_{\alpha, \beta} A(\alpha) S[\alpha|\beta] A(\beta)$.
  This sum is supposed to be the amplitude.
  If it's not symmetric, the kernel $S$ might be wrong.
  
  The Kernel I used: `klt_momentum_kernel_6pt`.
  Formula:
  $S[i_1 i_2 i_3 | j_1 j_2 j_3] = \prod (\dots)$.
  It depends on the choice of pivot "1" (index 0).
  Maybe I need to symmetrize over pivots?
  Or maybe the formula is just for specific helicities?
  Gravity MHV is helicity independent (except overall factor).
  
  **New Plan:**
  1. Analyze $\Delta N$.
  2. Try to find the "missing term".
  
  Let's run `analyze_numerator_diff.sage`.





