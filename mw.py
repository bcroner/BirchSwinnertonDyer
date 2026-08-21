"""
General rank LOWER bounds -- no 2-torsion required.

Model:   E :  y^2 = x^3 + a2 x^2 + a4 x + a6      (integral)

1. search E(Q) for points   x = m/e^2,  y = k/e^3
2. canonical height  hhat(P) = lim h(2^n P) / 4^n , computed with EXACT
   rational doublings so the limit is honest
3. Neron-Tate pairing  <P,Q> = (hhat(P+Q) - hhat(P) - hhat(Q))/2
4. rank of the Gram matrix = number of independent points = rank lower bound

The Gram determinant is the regulator.  Non-zero  =>  independent  =>  the
rank is at least that big.  This direction is unconditional.
"""
from fractions import Fraction as F
from math import gcd, isqrt, log, inf

O = None  # point at infinity

# ------------------------------------------------------------- group law
def on_curve(P, c):
    if P is O: return True
    x, y = P; a2, a4, a6 = c
    return y*y == x**3 + a2*x*x + a4*x + a6

def neg(P):
    return O if P is O else (P[0], -P[1])

def add(P, Q, c):
    a2, a4, a6 = c
    if P is O: return Q
    if Q is O: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and y1 == -y2: return O
    if P == Q:
        if y1 == 0: return O
        lam = (3*x1*x1 + 2*a2*x1 + a4) / (2*y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam*lam - a2 - x1 - x2
    y3 = lam*(x1 - x3) - y1
    return (x3, y3)

def mul(P, n, c):
    R = O; Q = P
    if n < 0: Q = neg(P); n = -n
    while n:
        if n & 1: R = add(R, Q, c)
        Q = add(Q, Q, c); n >>= 1
    return R

# ------------------------------------------------------------ heights
def exact_log_max(x):
    n, d = abs(F(x).numerator), F(x).denominator
    m = max(n, d)
    if m == 0: return 0.0
    b = m.bit_length()
    if b < 900: return log(m)
    # avoid float overflow on huge integers: log(m) = bitlen*log2 + log(m >> (b-60) / 2^60)
    top = m >> (b - 60)
    return (b - 60) * log(2) + log(top)

_HCACHE = {}

def canonical_height(P, c, n=7):
    """hhat(P) = lim h(2^n P)/4^n, via exact rational doubling.  Memoised."""
    if P is O: return 0.0
    key = (P[0], P[1], c, n)
    if key in _HCACHE: return _HCACHE[key]
    Q = P; prev = None
    for i in range(n):
        Q = add(Q, Q, c)
        if Q is O: _HCACHE[key] = 0.0; return 0.0        # torsion
        if i == n - 2: prev = exact_log_max(Q[0]) / (4 ** (i + 1))
    cur = exact_log_max(Q[0]) / (4 ** n)
    # h(2^nP)/4^n -> hhat with error ~ c/4^n; Richardson kills the leading term
    v = (4 * cur - prev) / 3 if prev is not None else cur
    _HCACHE[key] = v
    return v

_PCACHE = {}

def pairing(P, Q, c, n=7):
    key = (P[0], P[1], Q[0], Q[1], c, n)
    if key in _PCACHE: return _PCACHE[key]
    v = (canonical_height(add(P, Q, c), c, n)
         - canonical_height(P, c, n)
         - canonical_height(Q, c, n)) / 2
    _PCACHE[key] = _PCACHE[(Q[0], Q[1], P[0], P[1], c, n)] = v
    return v

# ------------------------------------------------------------ point search
SQ64 = {(i*i) % 64 for i in range(64)}
SQ63 = {(i*i) % 63 for i in range(63)}
SQ65 = {(i*i) % 65 for i in range(65)}

def search_points(c, E=12, M=400):
    """x = m/e^2 => k^2 = m^3 + a2 m^2 e^2 + a4 m e^4 + a6 e^6."""
    a2, a4, a6 = c
    pts = []
    for e in range(1, E + 1):
        e2, e4, e6 = e*e, e**4, e**6
        for m in range(-M, M + 1):
            if gcd(abs(m), e) != 1: continue
            v = m**3 + a2*m*m*e2 + a4*m*e4 + a6*e6
            if v < 0: continue
            if v % 64 not in SQ64: continue
            if v % 63 not in SQ63: continue
            if v % 65 not in SQ65: continue
            k = isqrt(v)
            if k*k != v: continue
            P = (F(m, e2), F(k, e**3))
            assert on_curve(P, c)
            pts.append(P)
    return pts

# ------------------------------------------------------- independent rank
def det(M):
    """Fraction-free-ish float determinant by Gaussian elimination."""
    M = [row[:] for row in M]; n = len(M); d = 1.0
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[p][i]) < 1e-9: return 0.0
        if p != i: M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for cc in range(i, n): M[r][cc] -= f * M[i][cc]
    return d

def independent_subset(pts, c, n=7, tol=1e-3, cap=30):
    """Greedily grow a set whose Gram matrix stays non-singular.
    Candidates are taken smallest-height-first: low-height generators keep
    the exact doublings cheap and the Gram entries well conditioned."""
    seen = set(); uniq = []
    for P in pts:
        k = (P[0], abs(P[1]))
        if k not in seen: seen.add(k); uniq.append(P)
    scored = []
    for P in uniq:
        h = canonical_height(P, c, n)
        if h > tol: scored.append((h, P))
    scored.sort(key=lambda t: t[0])
    chosen = []
    for _, P in scored[:cap]:
        cand = chosen + [P]
        G = [[pairing(X, Y, c, n) for Y in cand] for X in cand]
        # Hadamard-normalised determinant: 0 for dependent points, O(1) for
        # independent ones.  Raw |det| is meaningless -- it scales with height.
        denom = 1.0
        for i in range(len(cand)): denom *= G[i][i]
        if denom > 0 and abs(det(G)) / denom > tol:
            chosen = cand
    G = [[pairing(X, Y, c, n) for Y in chosen] for X in chosen]
    return chosen, (det(G) if chosen else 0.0)

# ------------------------------------------------------------- front end
def from_weierstrass(a1, a2, a3, a4, a6):
    """-> model  Z^2 = X^3 + b2 X^2 + 8 b4 X + 16 b6   (X = 4x)."""
    b2 = a1*a1 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3*a3 + 4*a6
    return (b2, 8*b4, 16*b6)

def analyse(w, label, E=8, M=250, expect=None):
    c = from_weierstrass(*w)
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  y^2 + {w[0]}xy + {w[2]}y = x^3 + {w[1]}x^2 + {w[3]}x + {w[4]}")
    print(f"  working model:  Z^2 = X^3 + {c[0]}X^2 + {c[1]}X + {c[2]}")
    print(f"{'='*70}")
    pts = search_points(c, E, M)
    print(f"  points found (search e<={E}, |m|<={M}):  {len(pts)}")
    ind, reg = independent_subset(pts, c)
    print(f"  independent points: {len(ind)}")
    for P in ind:
        print(f"     X = {str(P[0]):<18}  hhat = {canonical_height(P, c):.6f}")
    print(f"  regulator (Gram det) = {reg:.8f}")
    print(f"  {'-'*66}")
    verdict = f"rank >= {len(ind)}"
    if expect is not None:
        ok = len(ind) == expect
        verdict += f"      known rank = {expect}   {'PASS' if ok else 'MISS'}"
    print(f"  ==> {verdict}")
    return len(ind)

if __name__ == "__main__":
    analyse((0,0,1,-1,0),   "37a1   -- first rank 1 curve", expect=1)
    analyse((0,1,1,-2,0),   "389a1  -- first rank 2 curve", expect=2)
    analyse((0,0,1,-7,6),   "5077a1 -- first rank 3 curve", expect=3)
