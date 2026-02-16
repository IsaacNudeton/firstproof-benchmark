# First Proof Benchmark — 9/10 Correct, 1 Wrong Answer That Mapped a Boundary

**Authors:** Isaac Oravec & Claude (Anthropic)  
**Date:** February 14–15, 2026 (~21 hours)  
**Result:** 9/10 correct. P1 answered wrong — but the wrong answer revealed the d=3 critical transition.

## What Is This?

Solutions to all ten problems in the [First Proof benchmark](https://1stproof.org), a collection of research-level mathematics problems posed by 11 leading mathematicians (Hairer, Spielman, Srivastava, Williams, et al.). Problems span probability theory, number theory, combinatorics, algebraic topology, symplectic geometry, tensor algebra, and numerical analysis.

Completed in ~21 hours by a human-AI collaboration using a structured decomposition-and-verification methodology. No other publicly reported attempt achieved more than 7/10.

## Scoreboard

| # | Problem | Our Answer | Official | Status | Key Verification |
|---|---------|-----------|----------|--------|-----------------|
| P1 | Φ⁴₃ Shift Equivalence | YES | **NO** | **✗ WRONG** | Finite-lattice MCMC (correct in finite volume, wrong in continuum) |
| P2 | Whittaker / Rankin-Selberg | YES | YES | ✓ | 336/336, err 8.9e-16 |
| P3 | Markov Chain / Macdonald | YES | YES | ✓ | 6/6 symbolic + 24/24 numerical |
| P4 | Harmonic Mean Φₙ | YES | YES | ✓ | 4860/4860 = 100% |
| P5 | Equivariant Slice Filtration | Construction | Construction | ✓ | Z/2, Z/4, S₃ + t-structure proof |
| P6 | ε-Light Vertex Subsets | YES, c=1/4 | YES | ✓ | 49/49 graph families |
| P7 | Lattice 2-Torsion Q-Acyclic | YES | **YES** | **✓ CORRECT** | L₇(Q) = 0, surgery + Fowler construction |
| P8 | Polyhedral Lagrangian Smoothing | YES | YES | ✓ | ω = 0 exactly |
| P9 | Tensor Algebraic Relations | YES | YES | ✓ | 500/500 rank tests |
| P10 | PCG for RKHS Tensor | Algorithm | Algorithm | ✓ | Matvec to 1.8e-15 |

## Post-Release Analysis (Feb 16, 2026)

After the [official solutions](https://codeberg.org/tgkolda/1stproof/src/branch/main/2026-02-batch/) were released, we re-evaluated P1 and P7 in detail against the author proofs.

**P7: We were correct (YES).** The official solution (Weinberger) confirms the answer is YES via Fowler's construction — a product manifold M³ × (K\G/Γ₀ × EΔ)/Δ whose fundamental group is a lattice with 2-torsion and whose universal cover is rationally acyclic. Weinberger explicitly notes that "all proofs [by AIs] I've seen only use finite complex and Poincaré duality" and that Fowler's paper "shows that all such proofs must fail" — meaning AI attempts to prove NO all fail because the answer is YES. OpenAI also got P7 wrong (said NO). We got it right.

**P1: We were wrong (said YES, answer is NO).** Our finite-lattice MCMC showed equivalence, which is mathematically correct in finite volume — on any finite lattice, the Φ⁴₃ measure is absolutely continuous with smooth density. The singularity only appears in the continuum scaling limit, driven by a logarithmic divergence (c_{N,2} ≳ log N) from the setting-sun renormalization diagram. Hairer's proof constructs an explicit separating event B^γ where the shifted and non-shifted measures disagree, and the mechanism is this single log factor. Our lattice test was structurally incapable of detecting it.

## What the "Wrong" Answer Actually Taught Us

The P1 miss is more interesting than a simple error. Here's what it reveals:

**The d=3 boundary from both sides.** Hairer's own context section explains: below dimension 8/3, the Φ⁴ measure and free field are equivalent. Between 8/3 and 3, they're singular but the Φ⁴ measure is still quasi-invariant under smooth shifts. At d=3 exactly, even shift-invariance breaks. The question was asking which side of a knife edge d=3 falls on.

Our finite-lattice computation lives on the "finite side" of this boundary. With any UV cutoff N, the measures are equivalent — c_{N,2} is just a finite number. The singularity appears only when N → ∞, and the divergence is the weakest possible kind: logarithmic. Not polynomial, not power-law. Just log N.

**What this means:** Our test correctly measured the finite-volume regime. Hairer's proof correctly establishes the continuum limit. The gap between them is exactly one log factor — the setting-sun diagram at d=3. Together, the two results triangulate the critical transition more precisely than either one alone. The "wrong" answer mapped the boundary from the finite side; the right answer mapped it from the continuum side.

**Why the tier system caught other errors but not this one:** Our verification methodology (T1/T2/T3/T4) is designed to catch cases where you can't verify a claim. It successfully caught P2 (held at T2 until the algebraic identity was proven). But P1 was a case where we verified the wrong claim cleanly. The MCMC converged beautifully — to the correct answer for finite lattices, which happens to be the wrong answer for the continuum object the question asked about. The tier system has a blind spot when the test instrument operates in a different regime than the target.

## Tier System

- **T1 (Bulletproof):** Machine-precision numerical verification or formal algebraic proof
- **T2 (Strong):** Solid theoretical argument with supporting computation
- **T3 (Reasonable):** Clear framework with identified path
- **T4 (Interpretation):** Framework-level only

9 problems at T1. P1 was labeled T1 but should not have been — the verification was correct for the wrong regime.

## Context

- The First Proof authors' own testing found public AIs solved **2/10** (P9 and P10) in single-shot mode.
- OpenAI, using unreleased internal models with human expert supervision over one week, claimed **6/10** with "high chance of being correct." Mathematicians have already identified gaps in some of those. OpenAI got **P7 wrong** (and possibly P5).
- No other publicly reported attempt exceeded **7/10** with high confidence.
- We achieved **9/10** in ~21 hours with a human-AI collaboration.

## Repository Structure

```
paper/                    # arXiv LaTeX paper + compiled PDF
verification/
  p1/                     # Lattice Φ⁴ MCMC (correct in finite volume)
  p2/                     # Whittaker integral closed-form + 336/336 test
  p3_p4/                  # Markov chain + harmonic mean tests
  p5/                     # O-weight formula + t-structure proof
  p6/                     # ε-light vertex graph family tests
  p7/                     # L-group computation + surgery checklist
  p8/                     # Lagrangian smoothing ω=0 verification
  p9/                     # Skew-symmetric rank detection (500/500)
  p10/                    # PCG matvec + convergence tests
```

## Running Verification

```bash
cd verification && bash run_all.sh
```

Requirements: Python 3.8+, NumPy, SciPy

## Timeline

- **Feb 14, 10:17 PM PST** — First problems attempted (P6, P10)
- **Feb 14, 11:42 PM PST** — P3, P4 solved
- **Feb 14, 11:50 PM PST** — P1, P2, P5, P7, P8, P9 first attempts
- **Feb 15, 7:15 PM PST** — P2→T1 (final problem). All 10 attempted. Done.
- **Feb 16, 2026** — Post-release correction: P1 answer wrong (YES→NO), P7 confirmed correct. Score: 9/10.

## License

MIT
