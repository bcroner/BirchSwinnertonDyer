"""
ks2.py -- K(S,2) = { a in K*/K*^2 : v_q(a) even for every prime q not in S }.

Ideal factorisation is avoided by a content trick.  Write alpha = c * alpha'
where c is the largest RATIONAL integer dividing alpha.  Then alpha' is not
divisible by any rational prime, so for a split prime q = QQ' it can lie in
at most one of Q, Q'.  Hence, with j = v_q(c) and m = v_q(N(alpha)) - 2j:

    q inert     : v_Q(alpha) = j            -> need j even
    q ramified  : (q) = Q^2, v_Q = 2j + m   -> need m even
    q split     : {v_Q, v_Q'} = {j+m, j}    -> need j and m both even

All three are parity checks on ordinary integer valuations.
"""
from math import gcd, isqrt
from quadfield import Quad, squarefree_part
from descent import factor

def vq(n, q):
    if n == 0: return 10**9
    v = 0
    while n % q == 0: n //= q; v += 1
    return v

def splitting(K, q):
    """'ramified', 'split' or 'inert' for the rational prime q in K."""
    D = K.disc
    if D % q == 0: return "ramified"
    if q == 2:
        return "split" if D % 8 == 1 else "inert"
    return "split" if pow(D % q, (q-1)//2, q) == 1 else "inert"

def n_primes_above(K, S):
    n = 0
    for q in S:
        t = splitting(K, q)
        n += 2 if t == "split" else 1
    return n

def in_KS2(K, a, S):
    u, v = a
    N = K.norm(a)
    if N == 0: return False
    c = gcd(abs(u), abs(v))
    qs = set(factor(abs(N))) | (set(factor(c)) if c else set())
    for q in qs:
        if q in S: continue
        j = vq(c, q) if c else 0
        m = vq(abs(N), q) - 2*j
        t = splitting(K, q)
        if t == "inert":
            if j % 2: return False
        elif t == "ramified":
            if m % 2: return False
        else:
            if j % 2 or m % 2: return False
    return True

def same_class(K, a, b):
    """a ~ b in K*/K*^2  <=>  a*b is a square in K*."""
    return K.is_square_in_K(K.mul(a, b))

def compute_KS2(K, S, verbose=False):
    """Enumerate K(S,2) for IMAGINARY quadratic K.

    Every class has a representative alpha with (alpha) = (S-part) * b^2 and
    b chosen of norm <= Minkowski bound, so
        |N(alpha)| <= prod_{Q in S_K} N(Q)  *  M_K^2
    which is a PROVEN bound, not a guessed box."""
    assert K.d < 0, "imaginary quadratic only (real needs fundamental units)"
    Snorm = 1
    for q in S:
        t = splitting(K, q)
        Snorm *= (q*q if t == "inert" else q) if t != "split" else q*q
    M = K.minkowski()
    B = int(Snorm * M * M) + 4
    reps = []
    for a in K.elements_of_norm_upto(B):
        if not in_KS2(K, a, S): continue
        if any(same_class(K, a, r) for r in reps): continue
        reps.append(a)
    return reps, B

def is_group(K, reps):
    """Sanity: the representatives must close under multiplication mod squares."""
    for a in reps:
        for b in reps:
            p = K.mul(a, b)
            if not any(same_class(K, p, r) for r in reps): return False
    return True
