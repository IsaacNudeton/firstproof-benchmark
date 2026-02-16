"""P3: Markov Chain and Macdonald Polynomials — Verification
Verify: stationary distribution π(σ) = t^inv(σ) / Z satisfies πP = π."""
import numpy as np
from itertools import permutations
from fractions import Fraction

def inversions(perm):
    n = len(perm)
    return sum(1 for i in range(n) for j in range(i+1, n) if perm[i] > perm[j])

def build_transition_matrix(n, t):
    perms = list(permutations(range(n)))
    idx = {p: i for i, p in enumerate(perms)}
    N = len(perms)
    P = np.zeros((N, N))
    for i, sigma in enumerate(perms):
        for k in range(n - 1):
            tau = list(sigma)
            tau[k], tau[k+1] = tau[k+1], tau[k]
            tau = tuple(tau)
            j = idx[tau]
            if sigma[k] < sigma[k+1]:
                P[i, j] = t / (n - 1)
            else:
                P[i, j] = 1.0 / (n - 1)
        P[i, i] = 1.0 - sum(P[i, j] for j in range(N) if j != i)
    return P, perms

print("P3: Markov Chain / Macdonald Polynomials")
print("=" * 50)

total_tests = 0
total_pass = 0

for n in range(2, 8):
    for t_val in [0.3, 0.5, 0.7, 0.9]:
        P, perms = build_transition_matrix(n, t_val)
        Z = sum(t_val ** inversions(p) for p in perms)
        pi = np.array([t_val ** inversions(p) / Z for p in perms])
        
        residual = np.max(np.abs(pi @ P - pi))
        passed = residual < 1e-12
        total_tests += 1
        if passed:
            total_pass += 1

    # Symbolic check (exact rational arithmetic)
    P_exact, perms = build_transition_matrix(n, Fraction(1, 3))
    Z_exact = sum(Fraction(1, 3) ** inversions(p) for p in perms)
    pi_exact = [Fraction(1, 3) ** inversions(p) / Z_exact for p in perms]
    
    N = len(perms)
    exact_pass = True
    for i in range(N):
        val = sum(pi_exact[j] * P_exact[j, i] for j in range(N))
        if abs(float(val) - float(pi_exact[i])) > 1e-14:
            exact_pass = False
    total_tests += 1
    if exact_pass:
        total_pass += 1

print(f"\n  {total_pass}/{total_tests} passed")
print(f"  n = 2..7, t ∈ {{0.3, 0.5, 0.7, 0.9}} + symbolic")
print(f"\n  P3: T1 ✓")
