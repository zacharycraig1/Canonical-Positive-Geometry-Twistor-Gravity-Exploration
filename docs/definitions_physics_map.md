# Physics Definitions & Map Standards

This document serves as the authoritative source for spinor conventions, variable definitions, and the target gravity formula used in this repository.

## 1. Spinor Conventions

We use standard spinor helicity formalism variables:
- Holomorphic spinors: $\lambda_i^\alpha \in \mathbb{C}^2$
- Anti-holomorphic spinors: $\tilde\lambda_i^{\dot\alpha} \in \mathbb{C}^2$
- Lorentz invariant brackets:
  - $\langle ij \rangle = \epsilon_{\alpha\beta} \lambda_i^\alpha \lambda_j^\beta = \lambda_i^1 \lambda_j^2 - \lambda_i^2 \lambda_j^1$
  - $[ij] = \epsilon_{\dot\alpha\dot\beta} \tilde\lambda_i^{\dot\alpha} \tilde\lambda_j^{\dot\beta} = \tilde\lambda_i^{\dot 1} \tilde\lambda_j^{\dot 2} - \tilde\lambda_i^{\dot 2} \tilde\lambda_j^{\dot 1}$

## 2. Reference Spinors and Weights

For the mapping to the forest polynomial, we introduce two auxiliary reference spinors $x$ and $y$ (distinct from the external particle spinors).

The **weight** for particle $i$ is defined as:
$$
C_i = \langle i x \rangle \langle i y \rangle
$$

## 3. Edge Variables ($z_{ij}$)

The edge variables $z_{ij}$ used in the forest polynomial $F_{n,R}(z)$ are defined for $1 \le i < j \le n$ as:

$$
z_{ij} = \frac{[ij]}{\langle ij \rangle} \, C_i \, C_j
$$

## 4. The Verified Identity (Target)

The target identity relating the MHV gravity amplitude to the forest polynomial is:

$$
M_n^{\mathrm{MHV}} = (-1)^{n-1} \, \langle x y \rangle^8 \; \frac{F_{n,R}(z)}{\mathcal{N}_R \; \prod_{k \notin R} C_k^2}
$$

**Notes on Normalization:**
- The factor $\langle xy \rangle^8$ (sometimes denoted $\langle ab \rangle^8$ in literature) carries the helicity weight.
- $F_{n,R}(z)$ is the generating polynomial of forests where each component contains exactly one root from the set $R$.
- $\mathcal{N}_R$ is a normalization factor depending on the root set (usually related to the determinant of the Laplacian minor).
- The term $\prod_{k \notin R} C_k^2$ accounts for the weights of non-root vertices.

This formula matches the structure of the **Matrix-Tree Theorem** for weighted Laplacians (Hodges 2012, Feng & He 2012), where the determinant (amplitude) is expressed as a sum over forests.



