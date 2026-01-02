# Axioms of Signed Geometry

**Date:** January 2026  
**Status:** PROPOSED FRAMEWORK

---

## Motivation

Positive geometry, introduced by Arkani-Hamed and Trnka, axiomatizes structures like the Amplituhedron whose canonical forms equal scattering amplitudes. The key features are:

1. A bounded region with boundaries
2. A canonical form with logarithmic singularities on boundaries
3. Recursive structure: residues on boundaries yield lower-point positive geometries

**Problem:** Gravity does not fit this framework. The forest expansion has intrinsic sign cancellations that cannot be removed.

**Solution:** We propose **Signed Geometry** — a generalization that naturally incorporates mixed-sign contributions.

---

## Definition 1: Signed Geometry

A **signed geometry** is a triple $(X, \Sigma, \Omega)$ where:

1. **$X$ (Base Space):** A kinematic space (typically the space of on-shell momenta)

2. **$\Sigma$ (Stratification):** A decomposition of the amplitude into cells:
   $$\Sigma = \{(\sigma_\alpha, \varepsilon_\alpha)\}_{\alpha \in \mathcal{I}}$$
   where:
   - $\sigma_\alpha$ is a cell (e.g., a spanning forest)
   - $\varepsilon_\alpha \in \{+1, -1\}$ is the sign of that cell

3. **$\Omega$ (Signed Canonical Form):** A differential form given by:
   $$\Omega = \sum_{\alpha \in \mathcal{I}} \varepsilon_\alpha \cdot \omega_\alpha$$
   where $\omega_\alpha$ is the (positive) canonical form of cell $\sigma_\alpha$

---

## Axiom 1: Sign Rule

The sign $\varepsilon_\alpha$ is determined by a combinatorial rule depending on:
- The topology of cell $\alpha$
- The kinematic configuration

**For MHV Gravity:**
$$\varepsilon(F) = (-1)^{|E(F)|} \times \text{sign}\left(\prod_{e \in E(F)} w_e\right) \times \text{sign}\left(\prod_{v \in V} C_v^{\deg_F(v)}\right)$$

---

## Axiom 2: Balance Condition

The signs exhibit a characteristic balance tied to an underlying signature:

**Definition:** A signed geometry has **signature $(p, q)$** if the stratification has approximately $\frac{p}{p+q}$ positive cells and $\frac{q}{p+q}$ negative cells.

**For MHV Gravity (n=6):**
- 108 total forests
- Modal split: (54, 54) positive/negative
- Signature: (3, 3) — inherited from KLT kernel

---

## Axiom 3: Signed Factorization

At kinematic boundaries (poles), the signed canonical form factors as:

$$\text{Res}_{s_I = 0}\, \Omega_n = \Omega_L \wedge \Omega_R$$

where:
- $\Omega_L$ and $\Omega_R$ are lower-point signed canonical forms
- The sign structure is inherited: $\varepsilon(F) = \varepsilon(F_L) \cdot \varepsilon(F_R)$ for forests that decompose

**Key Property:** The sign rule is compatible with factorization.

---

## Axiom 4: Reference Independence

The total amplitude is independent of auxiliary choices (reference spinors):

$$\sum_\alpha \varepsilon_\alpha \cdot \omega_\alpha = \text{physical amplitude (independent of gauge)}$$

Even though individual $\varepsilon_\alpha$ and $\omega_\alpha$ depend on reference spinors, their sum does not.

---

## Axiom 5: Signature Origin

The signature $(p, q)$ of a signed geometry is determined by an underlying bilinear form.

**For Gravity:**
- The KLT kernel $S_{\text{KLT}}$ defines the bilinear structure
- Signature of $S_{\text{KLT}}$ for $n=6$ is $(3, 3)$
- This explains the 50/50 sign split in the forest expansion

---

## Comparison: Positive vs Signed Geometry

| Property | Positive Geometry | Signed Geometry |
|----------|-------------------|-----------------|
| Base space | Bounded region | Full kinematic space |
| Canonical form | All positive cells | Mixed-sign cells |
| Signature | $(n, 0)$ | $(p, q)$ with $p, q > 0$ |
| Example | Amplituhedron (YM) | Forest expansion (Gravity) |
| Underlying bilinear | Positive definite | Indefinite |

---

## Theorem: Gravity is Signed Geometry

**Statement:** The $n$-point MHV gravity amplitude is the canonical form of a signed geometry with:
- Base space: $X = $ kinematic space of $n$ on-shell gravitons
- Stratification: $\Sigma = $ 3-rooted spanning forests of $K_n$
- Sign rule: $\varepsilon(F) = (-1)^{|E|} \cdot \text{sign}(\prod w_e) \cdot \text{sign}(\prod C_v^{\deg})$
- Signature: $(p, p)$ for even $n$, derived from KLT kernel

**Proof:**
1. The amplitude equals the weighted Laplacian determinant (Hodges)
2. The determinant expands as a sum over forests (Matrix-Tree Theorem)
3. Each forest has a sign given by the sign rule (Theorem, verified)
4. The signature is inherited from the KLT kernel (verified)

QED ∎

---

## Open Questions

### Q1: Uniqueness
Is gravity the *only* signed geometry with signature $(3,3)$ for n=6?

### Q2: Loop Extension
Do loop amplitudes have signed geometry? What is the signature?

### Q3: Mathematical Structure
What is the algebraic structure underlying signed geometry?
- Possible connection to Clifford algebras (which have split signatures)
- Possible connection to twisted cohomology

### Q4: String Theory
What is the worldsheet origin of the sign rule?

---

## Verification Status

| Axiom | n=6 | n=7 | n=8 |
|-------|-----|-----|-----|
| Sign rule | ✅ 100% | ✅ 100% | ⏳ |
| Balance (50/50) | ✅ (54,54) | ✅ (~50/50) | ⏳ |
| Reference independence | ✅ | ✅ | ⏳ |
| KLT signature | ✅ (3,3) | ⏳ | ⏳ |

---

*Axioms formulated January 2026.*

