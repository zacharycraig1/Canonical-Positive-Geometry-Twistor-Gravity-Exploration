# Axioms for Signed Geometry: A Mathematical Framework

## Overview

**Signed geometry** is a generalization of positive geometry appropriate for theories whose scattering amplitudes arise from bilinear forms with **indefinite signature**.

## Definition: Signed Geometry

A **signed geometry** is a tuple (X, Omega, epsilon) where:
1. X is a geometric space (polytope, variety, or stratified space)
2. Omega is a meromorphic differential form (the canonical form)
3. epsilon: Cells(X) -> {+1, -1} is a sign function

such that: Omega = Sum_sigma epsilon(sigma) * omega(sigma)

## Axiom 1: Sign Coherence
The sign function must be coherent with boundaries.

## Axiom 2: Factorization Compatibility
If X factorizes as X_L x X_R, then: epsilon(sigma) = epsilon(sigma_L) x epsilon(sigma_R)

## Axiom 3: Signature Preservation
The signature (p, q) = (|C_+|, |C_-|) is an invariant of the signed geometry.

## Comparison: Positive vs Signed Geometry

| Property | Positive Geometry | Signed Geometry |
|----------|-------------------|-----------------|
| Sign function | epsilon = +1 | epsilon in {+1, -1} |
| Signature | (N, 0) | (p, q) with p,q > 0 |
| Example | Yang-Mills | Gravity |
| KLT role | N/A | Determines signature |

## The Gravity Instance

For n-point MHV gravity:
- Space: X = kinematic space with forest polytope
- Form: Omega = Hodges determinant = Sum_F epsilon(F) omega(F)
- Signs: epsilon(F) = sign(Prod w_e) x sign(Prod C_v^deg)
- Signature: (54, 54) for n=6, matching KLT's (3,3)

## Open Questions

1. Is the sign function unique?
2. Which theories have signed vs positive geometry?
3. Does signed geometry persist at loop level?
