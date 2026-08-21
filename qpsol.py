"""
qpsol.py -- decide  z^2 = g(u,v)  over Q_p  for a binary quartic g, odd p.

The two-variable refinement branches p^2 per level, which is hopeless at
p = 11 or 23.  Reduce to one variable instead:

  (u,v) primitive  =>  v is a unit (scale v=1) or p|v and u is a unit (u=1).

Then decide  ?exists u in Z_p : mult * f(u) in (Q_p*)^2 u {0}  recursively:

  * strip the content: p^2 divides out freely (a square); a single p moves
    into `mult`, which therefore only ever takes the values 1 or p times a
    unit class.
  * for each residue r mod p, f(u) = f(r) mod p, so when v_p(f(r)) = 0 the
    square class is already decided for the whole residue class -- accept or
    discard it outright.
  * only when p | f(r) do we recurse, on f(r + p t).

Branching is p per level and the content division keeps the coefficients
bounded, so this terminates quickly where the 2-variable search cannot.
"""
from math import gcd

def _vp(n, p):
    if n == 0: return 10**9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

def _unit_sq(u, p):
    if p == 2: return u % 8 == 1        # unit square in Z_2 iff = 1 mod 8
    return pow(u % p, (p-1)//2, p) == 1

def _ev(f, x):
    r = 0
    for c in f: r = r*x + c
    return r

def _shift(f, r, p):
    """coefficients of f(r + p t), highest degree first."""
    n = len(f) - 1
    out = [0]*(n+1)
    # Horner-style synthetic expansion
    cur = list(f)
    for k in range(n+1):
        # evaluate and deflate at r
        acc = 0; newc = []
        for c in cur:
            newc.append(acc)
            acc = acc*r + c
        out[n-k] = acc
        cur = newc[1:] if len(newc) > 1 else []
        if not cur: break
    for k in range(n+1):
        out[n-k] *= p**k
    return out

def _step(p):
    """Residue modulus to branch on.  For odd p the square class of a unit is
    fixed by its residue mod p, so branching mod p suffices.  In Q_2 the class
    is fixed only mod 8, so we branch mod 8 and substitute u = r + 8t."""
    return 8 if p == 2 else p

def _sol1(f, mult, p, depth, maxdepth, seen):
    if depth > maxdepth: return None
    f = list(f)
    while f and f[0] == 0: f = f[1:]
    if not f: return True                      # f identically 0
    j = min(_vp(c, p) for c in f if c != 0)
    if j >= 2:
        f = [c // p**(2*(j//2)) for c in f]
        j = j - 2*(j//2)
    if j == 1:
        f = [c // p for c in f]
        mult = mult * p
    mv = _vp(mult, p)
    if mv >= 2: mult //= p**(2*(mv//2))
    key = (tuple(f), mult)
    if key in seen: return False
    seen = seen | {key}
    unknown = False
    st = _step(p)
    for r in range(st):
        val = _ev(f, r)
        if val == 0: return True               # z = 0 is a point
        m = _vp(val, p)
        if m == 0:
            tot = mult * val
            if _vp(tot, p) % 2 == 0 and _unit_sq(tot // p**_vp(tot, p), p):
                return True
            continue                            # class decided: discard
        res = _sol1(_shift(f, r, st), mult, p, depth+1, maxdepth, seen)
        if res is True: return True
        if res is None: unknown = True
    return None if unknown else False

def qp_soluble(g, p, maxdepth=40):
    """z^2 = g(u,v) over Q_p, g = (a,b,c,d,e) meaning a u^4 + ... + e v^4."""
    a, b, c, d, e = g
    # branch 1: v a unit -> f(u) = g(u,1)
    r1 = _sol1([a, b, c, d, e], 1, p, 0, maxdepth, frozenset())
    if r1 is True: return True
    # branch 2: p | v, u a unit -> g(1, p t)
    f2 = [e*p**4, d*p**3, c*p*p, b*p, a]      # g(1, p t)
    r2 = _sol1(f2, 1, p, 0, maxdepth, frozenset())
    if r2 is True: return True
    if r1 is None or r2 is None: return None
    return False
