"""
Problem 9 — Attack 3: The 2x2 Minor Approach

KEY INSIGHT I MISSED:
For fixed (γ,δ,k,l), define the matrix:
  M_{(α,i),(β,j)} = R^{αβγδ}_{ijkl} = λ_{αβγδ} · Q^{αβγδ}_{ijkl}

If λ is rank-1 in modes (1,2): λ_{αβγδ} = u_α · v_β · (w_γ x_δ)
Then M_{(α,i),(β,j)} = u_α · v_β · w_γ x_δ · Q^{αβγδ}_{ijkl}

The Q part: Q^{αβγδ}_{ijkl} = det[A^α(i,:); A^β(j,:); A^γ(k,:); A^δ(l,:)]
For fixed γ,δ,k,l, let c = A^γ(k,:), d = A^δ(l,:).
Q = det[a_{α,i}; b_{β,j}; c; d] where a_{α,i} = A^α(i,:), b_{β,j} = A^β(j,:)

This is BILINEAR in a_{α,i} and b_{β,j}. Specifically:
Q = Σ_{m,p} a_{α,i,m} · b_{β,j,p} · C_{mp}
where C_{mp} = cofactor involving c,d.

So Q^{αβγδ}_{ijkl} = a_{α,i}^T · C · b_{β,j}  (inner product through matrix C)

Now: M_{(α,i),(β,j)} = λ_{αβγδ} · (a_{α,i}^T · C · b_{β,j})

If λ = u⊗v⊗w⊗x:
M_{(α,i),(β,j)} = u_α v_β (w_γ x_δ) · (a_{α,i}^T · C · b_{β,j})
                = (u_α · a_{α,i})^T · (w_γ x_δ · C) · (v_β · b_{β,j})

Define: s_{α,i} = u_α · a_{α,i} (scaled row vector in R^4)
        t_{β,j} = v_β · b_{β,j}

Then M_{(α,i),(β,j)} = s_{α,i}^T · C' · t_{β,j}  where C' = (w_γ x_δ) · C

This is a BILINEAR form: M = S^T · C' · T
where S is 4 × 3n (columns are the s vectors)
and T is 4 × 3n (columns are the t vectors)

So rank(M) ≤ 4 (since C' is 4×4, or actually it might be smaller).

Wait — C is the cofactor matrix. For 4×4 determinant with rows (a,b,c,d):
det = Σ_{m} a_m · cofactor_m(b,c,d)
    = Σ_{m,p} a_m · b_p · cofactor_{mp}(c,d)

Actually det[a; b; c; d] as a function of a and b (with c,d fixed):
= Σ_σ sign(σ) a_{σ(1)} b_{σ(2)} c_{σ(3)} d_{σ(4)}

Let me think of it differently. The 4×4 determinant is the wedge product:
a ∧ b ∧ c ∧ d

As a bilinear function of a and b:
f(a,b) = a ∧ b ∧ c ∧ d

In coordinates: f(a,b) = Σ_{m<p} (a_m b_p - a_p b_m) · (c ∧ d)_{m̄p̄}
where (m̄,p̄) are the complement indices.

So f(a,b) = Σ_{m<p} (a_m b_p - a_p b_m) · C_{mp}

This is an ANTISYMMETRIC bilinear form in (a,b)!
f(a,b) = a^T · A_cd · b  where A_cd is antisymmetric (skew-symmetric).

A 4×4 skew-symmetric matrix has rank 0, 2, or 4.
For generic c,d it has rank 2 (since c∧d determines a 2-form, and 
the Hodge dual gives the complementary 2-plane).

So: rank(A_cd) = 2 for generic c,d.

Therefore: M = S^T · A_cd' · T has rank ≤ 2!

When λ is rank-1: the 3n × 3n matrix M has rank ≤ 2.
When λ is NOT rank-1: M can have higher rank.

THIS IS THE TEST!
"""
import numpy as np
np.random.seed(42)

n = 5

# Generate generic A matrices
A = [np.random.randn(3, 4) for _ in range(n)]

def det4(a, b, c, d):
    """4x4 determinant of row vectors a,b,c,d in R^4"""
    M = np.vstack([a, b, c, d])
    return np.linalg.det(M)

def build_M_matrix(A_list, n, gamma, delta, k, l, lam):
    """
    Build the 3n × 3n matrix M_{(α,i),(β,j)} = λ_{αβγδ} · Q^{αβγδ}_{ijkl}
    for fixed (γ,δ,k,l).
    """
    dim = 3 * n
    M = np.zeros((dim, dim))
    
    c = A_list[gamma][k, :]  # fixed row
    d = A_list[delta][l, :]  # fixed row
    
    for alpha in range(n):
        for i in range(3):
            row_idx = alpha * 3 + i
            a = A_list[alpha][i, :]
            for beta in range(n):
                for j in range(3):
                    col_idx = beta * 3 + j
                    b = A_list[beta][j, :]
                    
                    Q_val = det4(a, b, c, d)
                    
                    # λ is nonzero only when indices are not all identical
                    if alpha == beta == gamma == delta:
                        lam_val = 0
                    else:
                        lam_val = lam[alpha, beta, gamma, delta]
                    
                    M[row_idx, col_idx] = lam_val * Q_val
    
    return M

print("=" * 70)
print("P9 ATTACK 3: RANK TEST VIA SKEW-SYMMETRIC STRUCTURE")
print("=" * 70)

# Test 1: Rank-1 lambda
print("\n--- Rank-1 λ = u⊗v⊗w⊗x ---")
u = np.random.randn(n) * 2 + 0.5
v = np.random.randn(n) * 2 + 0.5
w = np.random.randn(n) * 2 + 0.5
x = np.random.randn(n) * 2 + 0.5
lam1 = np.einsum('a,b,c,d->abcd', u, v, w, x)

# Test multiple (γ,δ,k,l) choices
for gamma, delta, k, l in [(0,1,0,0), (1,2,1,2), (0,3,2,1), (2,4,0,2)]:
    if gamma == delta:
        continue
    M = build_M_matrix(A, n, gamma, delta, k, l, lam1)
    sv = np.linalg.svd(M, compute_uv=False)
    # Count significant singular values
    threshold = 1e-10 * sv[0] if sv[0] > 0 else 1e-10
    rank = np.sum(sv > threshold)
    print(f"  (γ={gamma},δ={delta},k={k},l={l}): rank(M) = {rank}, "
          f"sv[:5] = {sv[:5].round(6)}")

# Test 2: Random (non-rank-1) lambda
print("\n--- Random (non-rank-1) λ ---")
lam2 = np.random.randn(n, n, n, n)

for gamma, delta, k, l in [(0,1,0,0), (1,2,1,2), (0,3,2,1), (2,4,0,2)]:
    if gamma == delta:
        continue
    M = build_M_matrix(A, n, gamma, delta, k, l, lam2)
    sv = np.linalg.svd(M, compute_uv=False)
    threshold = 1e-10 * sv[0] if sv[0] > 0 else 1e-10
    rank = np.sum(sv > threshold)
    print(f"  (γ={gamma},δ={delta},k={k},l={l}): rank(M) = {rank}, "
          f"sv[:5] = {sv[:5].round(6)}")

# Test 3: Multiple rank-1 trials
print("\n--- Multiple rank-1 trials ---")
max_ranks = []
for trial in range(20):
    u = np.random.randn(n) * 3
    v = np.random.randn(n) * 3
    w = np.random.randn(n) * 3
    x = np.random.randn(n) * 3
    lam = np.einsum('a,b,c,d->abcd', u, v, w, x)
    
    max_rank = 0
    for gamma, delta, k, l in [(0,1,0,0), (1,2,1,2), (0,3,2,1)]:
        M = build_M_matrix(A, n, gamma, delta, k, l, lam)
        sv = np.linalg.svd(M, compute_uv=False)
        threshold = 1e-10 * sv[0] if sv[0] > 0 else 1e-10
        rank = np.sum(sv > threshold)
        max_rank = max(max_rank, rank)
    
    max_ranks.append(max_rank)
    status = "✓" if max_rank <= 2 else "✗"
    print(f"  Trial {trial+1}: max rank = {max_rank} {status}")

print(f"\n  All rank-1 trials have max rank ≤ 2: {all(r <= 2 for r in max_ranks)}")

# Test 4: Rank-2 lambda (should have higher rank)
print("\n--- Rank-2 λ ---")
u1, u2 = np.random.randn(n), np.random.randn(n)
v1, v2 = np.random.randn(n), np.random.randn(n)
w1, w2 = np.random.randn(n), np.random.randn(n)
x1, x2 = np.random.randn(n), np.random.randn(n)
lam_r2 = (np.einsum('a,b,c,d->abcd', u1, v1, w1, x1) + 
           np.einsum('a,b,c,d->abcd', u2, v2, w2, x2))

for gamma, delta, k, l in [(0,1,0,0), (1,2,1,2), (0,3,2,1), (2,4,0,2)]:
    if gamma == delta:
        continue
    M = build_M_matrix(A, n, gamma, delta, k, l, lam_r2)
    sv = np.linalg.svd(M, compute_uv=False)
    threshold = 1e-10 * sv[0] if sv[0] > 0 else 1e-10
    rank = np.sum(sv > threshold)
    print(f"  (γ={gamma},δ={delta},k={k},l={l}): rank(M) = {rank}, "
          f"sv[:5] = {sv[:5].round(6)}")

print("\n" + "=" * 70)
print("ANALYSIS")
print("=" * 70)
