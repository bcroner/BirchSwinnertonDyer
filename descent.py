"""
2-descent rank calculator for elliptic curves over Q.

Method: descent via 2-isogeny (Silverman, AEC X.4.9).

    E  : y^2 = x^3 + a x^2 + b x
    E' : Y^2 = X^3 - 2a X^2 + (a^2 - 4b) X

    rank E(Q) = log2|alpha(E(Q))| + log2|alpha'(E'(Q))| - 2

alpha sends (x,y) to the squarefree part of x.  For squarefree d | b, the
class d is in the image iff the homogeneous space

    N_d :  w^2 = d u^4 + a u^2 v^2 + (b/d) v^4

has a rational point.  Everywhere-locally-solvable N_d  ->  Selmer  ->  upper
bound on rank.  N_d with an actual rational point  ->  lower bound.
The difference between the two bounds is Sha[2].
"""
from math import gcd, isqrt

# --------------------------------------------------------------- arithmetic
def factor(n):
    n = abs(n); f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def squarefree_divisors(n):
    """All squarefree d (positive and negative) with d | n."""
    ps = sorted(factor(n))
    out = [1]
    for p in ps:
        out += [d * p for d in out]
    return sorted(set(out + [-d for d in out]), key=lambda z: (abs(z), z))

def squarefree_part(q):
    """Squarefree integer representing the rational q modulo squares."""
    num, den = q.numerator, q.denominator
    n = num * den
    s = 1
    for p, e in factor(n).items():
        if e % 2: s *= p
    return s if n > 0 else -s

def sqfree_int(n):
    """Squarefree integer representing n modulo squares."""
    if n == 0: return 0
    s = 1
    for p, e in factor(n).items():
        if e % 2: s *= p
    return s if n > 0 else -s

def saturate(classes):
    """Close a set of squarefree classes under multiplication mod squares.
    alpha's image and the Selmer group are both SUBGROUPS, so this is sound
    and it recovers classes the point search was too weak to witness."""
    G = set(classes) | {1}
    changed = True
    while changed:
        changed = False
        for x in list(G):
            for y in list(G):
                z = sqfree_int(x * y)
                if z not in G: G.add(z); changed = True
    return G

def vp(n, p):
    v = 0
    while n % p == 0: n //= p; v += 1
    return v, n

def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p-1)//2, p) == 1 else -1)

def is_padic_square(g, p):
    """Exact test: is the nonzero integer g a square in Q_p?"""
    m, u = vp(abs(g), p)
    if m % 2: return False
    if g < 0: u = -u
    return (u % 8 == 1) if p == 2 else legendre(u, p) == 1

def is_square(n):
    if n < 0: return False
    r = isqrt(n); return r * r == n

# ----------------------------------------------------- local solvability
def locally_solvable(d, a, e, p, node_budget=400_000, maxdepth=14):
    """
    Is  w^2 = d u^4 + a u^2 v^2 + e v^4  solvable over Q_p?

    Rigorous refinement: if (u,v) = (u0,v0) mod p^k then g(u,v) = g(u0,v0)
    mod p^k.  So once v_p(g(u0,v0)) = m < k, the valuation is pinned for
    every lift, and the unit part is known mod p^(k-m).  Knowing it mod p
    (p odd) or mod 8 (p=2) decides squareness outright.  Otherwise refine.
    """
    g = lambda u, v: d*u**4 + a*u*u*v*v + e*v**4
    need = 3 if p == 2 else 1
    unknown = False
    stack = [(u, v, 1) for u in range(p) for v in range(p) if (u or v)]
    seen = 0
    while stack:
        u0, v0, k = stack.pop()
        seen += 1
        if seen > node_budget: return None            # honest "don't know"
        val = g(u0, v0)
        if val == 0: return True                      # (u0,v0,0) is a point
        m, _ = vp(abs(val), p)
        if m < k and k - m >= need:
            if is_padic_square(val, p): return True
            continue                                  # this branch is dead
        if k >= maxdepth: unknown = True; continue
        pk = p ** k
        for s in range(p):
            for t in range(p):
                stack.append((u0 + s*pk, v0 + t*pk, k + 1))
    return None if unknown else False

NB = [400_000]
MD = [14]

def everywhere_locally_solvable(d, a, e, disc, verbose=False):
    """Real place plus every bad prime.  Good odd p >= 5 are automatic:
    good reduction + Hasse-Weil gives #C(F_p) >= p+1-2sqrt(p) > 0."""
    if not _real_solvable(d, a, e): return False
    bad = set(factor(2 * disc)) | {2}
    for p in sorted(bad):
        r = locally_solvable(d, a, e, p, node_budget=NB[0], maxdepth=MD[0])
        if r is False: return False
        if r is None: return None
    return True

def _real_solvable(d, a, e):
    """Does d u^4 + a u^2 v^2 + e v^4 take a non-negative value?"""
    if d > 0 or e > 0: return True
    if d == 0 or e == 0: return True
    # d,e < 0: need a > 0 and discriminant a^2 - 4de >= 0 for a positive bump
    return a > 0 and a * a - 4 * d * e >= 0

# ------------------------------------------------------------ point search
SQ64 = {(i*i) % 64 for i in range(64)}
SQ63 = {(i*i) % 63 for i in range(63)}
SQ65 = {(i*i) % 65 for i in range(65)}

def global_point(d, a, e, B=600):
    """Search for a rational point on w^2 = d u^4 + a u^2 v^2 + e v^4.
    Quadratic-residue sieve mod 64/63/65 kills ~99% of candidates before
    the expensive integer square root."""
    for v in range(0, B + 1):
        v2 = v * v; v4 = v2 * v2; ev4 = e * v4
        us = range(1, B + 1) if v else range(1, 2)
        for u in us:
            if v and gcd(u, v) != 1: continue
            u2 = u * u
            val = d*u2*u2 + a*u2*v2 + ev4
            if val < 0: continue
            if val % 64 not in SQ64: continue
            if val % 63 not in SQ63: continue
            if val % 65 not in SQ65: continue
            if is_square(val): return (u, v, isqrt(val))
    return None

# ------------------------------------------------------------ the descent
def alpha_image(a, b, B=600, verbose=False):
    """Returns (found_classes, selmer_classes, undetermined) for y^2=x^3+ax^2+bx."""
    disc = b * b * (a * a - 4 * b)
    found, selmer, undet = {1}, {1}, []
    for d in squarefree_divisors(b):
        if d == 1: continue
        if b % d: continue
        e = b // d
        pt = global_point(d, a, e, B)
        if pt: found.add(d); selmer.add(d); continue
        els = everywhere_locally_solvable(d, a, e, disc)
        if els is True: selmer.add(d)
        elif els is None: selmer.add(d); undet.append(d)
    found, selmer = saturate(found), saturate(selmer)
    assert found <= selmer, "image must sit inside Selmer"
    return found, selmer, undet

def log2size(s):
    """Exact log2 of the order of a subgroup of Q*/Q*^2."""
    n = len(s); k = n.bit_length() - 1
    assert 1 << k == n, f"not a 2-group: order {n} -- saturation failed"
    return k

def rank_bounds(a, b, B=600):
    """Full 2-isogeny descent.  Returns dict of results."""
    a2, b2 = -2 * a, a * a - 4 * b
    fE, sE, uE = alpha_image(a, b, B)
    fF, sF, uF = alpha_image(a2, b2, B)
    lo = log2size(fE) + log2size(fF) - 2
    hi = log2size(sE) + log2size(sF) - 2
    return dict(a=a, b=b, a2=a2, b2=b2,
                found_E=sorted(fE), selmer_E=sorted(sE), undet_E=uE,
                found_F=sorted(fF), selmer_F=sorted(sF), undet_F=uF,
                lo=max(lo, 0), hi=hi)
