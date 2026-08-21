"""
complete2.py -- COMPLETE 2-descent.  No search box, no reduction theory.

For E : y^2 = (x-e1)(x-e2)(x-e3) with all e_i rational, the descent map

    alpha : E(Q)/2E(Q) -> Q*/Q*^2 x Q*/Q*^2 ,  P |-> (x-e1, x-e2)

has image inside the subgroup supported on
    S = {-1} u {2} u {p : p | (e1-e2)(e1-e3)(e2-e3)} .

That subgroup is FINITE AND EXPLICIT -- it is just the squarefree integers
built from S -- so the candidate set is enumerated exactly, with no box and
no reduction theory.  This is why the etale algebra splitting matters: for
L = Q x Q x Q the S-unit group is elementary.

For each candidate (d1,d2), put d3 = squarefree(d1*d2).  The homogeneous
space is the intersection of two quadrics

    d1 u1^2 - d2 u2^2 = (e2-e1) w^2
    d1 u1^2 - d3 u3^2 = (e3-e1) w^2

so, writing A = d1 u1^2 - (e2-e1) w^2 and B = d1 u1^2 - (e3-e1) w^2, we need
some primitive (u1,w) with BOTH d2*A and d3*B squares.

    everywhere locally solvable  ->  the class lies in Sel_2
    a rational point             ->  the class lies in alpha(E(Q))

    dim Sel_2 = log2 #ELS ,   rank <= dim Sel_2 - 2   (since E(Q)[2]=(Z/2)^2)
"""
from math import isqrt, gcd
from descent import factor, is_padic_square, vp, sqfree_int

def support(e1, e2, e3):
    d = abs((e1-e2)*(e1-e3)*(e2-e3))
    return sorted(set(factor(d)) | {2})

def sqfree_classes(S):
    """All squarefree integers (both signs) built from the primes in S."""
    out = [1]
    for p in S: out += [d*p for d in out]
    return sorted(set(out + [-d for d in out]), key=lambda z: (abs(z), z))

# --------------------------------------------------------------- local test
def _pair_sq(d2, A, d3, B, p):
    """Are d2*A and d3*B both squares in Q_p?  (A,B nonzero integers.)"""
    if A == 0 and B == 0: return True
    if A == 0: return is_padic_square(d3*B, p)
    if B == 0: return is_padic_square(d2*A, p)
    return is_padic_square(d2*A, p) and is_padic_square(d3*B, p)

def p_solvable(d1, d2, d3, e1, e2, e3, p, maxdepth=22, budget=4_000_000):
    """Exhaustive p-adic test by valuation pinning: once v_p(F(u1,w)) = m < k,
    the valuation is fixed for every lift and the unit part is known mod
    p^(k-m); knowing it mod p (p odd) or mod 8 (p=2) decides squareness."""
    E2, E3 = e2-e1, e3-e1
    need = 3 if p == 2 else 1
    unknown = False
    stack = [(u, w, 1) for u in range(p) for w in range(p) if (u or w)]
    seen = 0
    while stack:
        u1, w, k = stack.pop()
        seen += 1
        if seen > budget: return None
        A = d1*u1*u1 - E2*w*w
        B = d1*u1*u1 - E3*w*w
        dead = False; all_decided = True
        for (dd, V) in ((d2, A), (d3, B)):
            if V == 0: continue          # exact zero: 0 is a square, condition met
            val = dd*V
            m, _ = vp(abs(val), p)
            if m < k and k - m >= need:
                if not is_padic_square(val, p):
                    dead = True; break   # PRUNE: this branch cannot ever work
            else:
                all_decided = False
        if dead: continue
        if all_decided: return True
        if k >= maxdepth: unknown = True; continue
        pk = p**k
        for s in range(p):
            for t in range(p):
                stack.append((u1 + s*pk, w + t*pk, k+1))
    return None if unknown else False

def real_solvable(d1, d2, d3, e1, e2, e3):
    E2, E3 = e2-e1, e3-e1
    for num in range(-400, 401):
        for den in (1, 2, 3, 5, 7):
            u1 = num/den
            A = d1*u1*u1 - E2
            B = d1*u1*u1 - E3
            if d2*A >= 0 and d3*B >= 0: return True
    # w = 0 branch: need d1*d2 > 0 and d1*d3 > 0
    return d1*d2 > 0 and d1*d3 > 0

def els(d1, d2, d3, e1, e2, e3, S):
    if not real_solvable(d1, d2, d3, e1, e2, e3): return False
    for p in sorted(set(S) | {2}):
        r = p_solvable(d1, d2, d3, e1, e2, e3, p)
        if r is False: return False
        if r is None: return None
    return True

_SQ64 = {(i*i) % 64 for i in range(64)}
_SQ63 = {(i*i) % 63 for i in range(63)}
_SQ65 = {(i*i) % 65 for i in range(65)}

def rational_point(d1, d2, d3, e1, e2, e3, B=1200):
    """Search the homogeneous space for a rational point.

    A projective point scales so that (u1,u2,u3,w) are coprime AS A QUADRUPLE.
    u1 and w may therefore share a factor -- requiring gcd(u1,w)=1 (an earlier
    bug) silently discarded valid points and cost half the image on n=37.
    """
    E2, E3 = e2-e1, e3-e1
    if d1*d2 > 0 and d1*d3 > 0:
        if isqrt(d1*d2)**2 == d1*d2 and isqrt(d1*d3)**2 == d1*d3:
            return ("infinity", 1, 0)
    for w in range(0, B+1):
        w2 = w*w; c2 = E2*w2; c3 = E3*w2
        us = range(-B, B+1) if w else range(1, 2)
        for u1 in us:
            t = d1*u1*u1
            x = d2*(t - c2)
            if x < 0: continue
            if x % 64 not in _SQ64: continue
            y = d3*(t - c3)
            if y < 0: continue
            if y % 64 not in _SQ64: continue
            if x % 63 not in _SQ63 or y % 63 not in _SQ63: continue
            if x % 65 not in _SQ65 or y % 65 not in _SQ65: continue
            if isqrt(x)**2 != x or isqrt(y)**2 != y: continue
            return (u1, w, t)
    return None

def saturate_pairs(pairs):
    """alpha(E(Q)) is a SUBGROUP of (Q*/Q*^2)^2 under componentwise
    multiplication mod squares, so closing the witnessed set under that
    operation is sound -- and it recovers classes the point search was too
    weak to exhibit.  (Same argument as the saturation in descent.py.)"""
    G = set(pairs) | {(1, 1)}
    changed = True
    while changed:
        changed = False
        for (a1, a2) in list(G):
            for (b1, b2) in list(G):
                z = (sqfree_int(a1*b1), sqfree_int(a2*b2))
                if z not in G: G.add(z); changed = True
    return G

def complete_descent(e1, e2, e3, verbose=False):
    S = support(e1, e2, e3)
    cands = sqfree_classes(S)
    sel, img, unk = [], [], []
    for d1 in cands:
        for d2 in cands:
            d3 = sqfree_int(d1*d2)
            r = els(d1, d2, d3, e1, e2, e3, S)
            if r is None:
                # SOUNDNESS: an undecided class must be counted INSIDE Selmer.
                # Dropping it shrinks |Sel_2| and hence the rank UPPER bound,
                # which is the unsafe direction -- it can report a rank that is
                # too small.  Counting it in can only make the bound weaker.
                unk.append((d1, d2)); sel.append((d1, d2)); continue
            if not r: continue
            sel.append((d1, d2))
            if rational_point(d1, d2, d3, e1, e2, e3): img.append((d1, d2))
    img = sorted(saturate_pairs(img))
    selset = set(sel)
    assert set(img) <= selset, "image must lie inside Selmer -- bug"
    return S, cands, sel, img, unk
