"""
Stage 3b: proper GL2(Z) equivalence of binary quartics, plus the
everywhere-locally-solvable filter.
"""
from math import isqrt, gcd
from quartic import I_inv, J_inv, disc
from descent import is_padic_square, vp, factor

def act(g, M):
    """g(u,v) -> g(al u + be v, ga u + de v).  Coefficient index = power of v."""
    (al, be), (ga, de) = M
    def pmul(P, Q):
        R = [0]*(len(P)+len(Q)-1)
        for i, x in enumerate(P):
            for j, y in enumerate(Q): R[i+j] += x*y
        return R
    tot = [0]*5
    for k, coef in enumerate(g):                 # term u^(4-k) v^k
        P = [1]
        for _ in range(4-k): P = pmul(P, [al, be])
        for _ in range(k):   P = pmul(P, [ga, de])
        for j, x in enumerate(P): tot[j] += coef*x
    return tuple(tot)

GENS = [((1,1),(0,1)), ((1,-1),(0,1)), ((1,0),(1,1)), ((1,0),(-1,1)),
        ((0,1),(1,0)), ((-1,0),(0,1)), ((1,0),(0,-1))]

def size(g): return sum(x*x for x in g)

def canonical(g, beam=40, rounds=30):
    """Beam-search the GL2(Z) orbit for the smallest representative, then
    take the lexicographic minimum among all forms of that minimal size."""
    frontier = {tuple(g)}
    best = {tuple(g)}
    bs = size(g)
    for _ in range(rounds):
        nxt = set()
        for h in frontier:
            for M in GENS:
                k = act(h, M)
                nxt.add(k)
        pool = sorted(nxt | best, key=size)[:beam]
        nb = size(pool[0])
        if nb < bs:
            bs = nb; best = set(pool); frontier = set(pool)
        else:
            best |= {p for p in pool if size(p) == bs}
            frontier = set(pool)
            if _ > 6: break
    minimal = [h for h in best if size(h) == bs]
    return min(minimal)

# ------------------------------------------------ everywhere locally solvable
def q_eval(g, u, v):
    a, b, c, d, e = g
    return a*u**4 + b*u**3*v + c*u*u*v*v + d*u*v**3 + e*v**4

def real_solvable(g):
    """Does g(u,v) take a value >= 0?"""
    for u, v in [(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(3,1),(1,3),(-1,2),(2,-1),
                 (5,1),(1,5),(3,2),(2,3),(7,2),(2,7),(4,1),(1,4)]:
        if q_eval(g, u, v) >= 0: return True
    for k in range(-4000, 4001):
        if q_eval(g, k, 1) >= 0: return True
    return False

def p_solvable(g, p, maxdepth=12, budget=300_000):
    """w^2 = g(u,v) over Q_p, by the same rigorous refinement as descent.py:
    once v_p(g(u0,v0)) = m < k the valuation is pinned for all lifts."""
    need = 3 if p == 2 else 1
    unknown = False
    stack = [(u, v, 1) for u in range(p) for v in range(p) if (u or v)]
    seen = 0
    while stack:
        u0, v0, k = stack.pop()
        seen += 1
        if seen > budget: return None
        val = q_eval(g, u0, v0)
        if val == 0: return True
        m, _ = vp(abs(val), p)
        if m < k and k - m >= need:
            if is_padic_square(val, p): return True
            continue
        if k >= maxdepth: unknown = True; continue
        pk = p**k
        for s in range(p):
            for t in range(p):
                stack.append((u0 + s*pk, v0 + t*pk, k+1))
    return None if unknown else False

def els(g, I, J):
    """Everywhere local solvability of z^2 = g(u,v).
    Uses the one-variable p-adic test (qpsol), which branches p per level
    instead of p^2 and divides out content, so it decides cases the
    two-variable refinement could only time out on."""
    from qpsol import qp_soluble
    if not real_solvable(g): return False
    D = 4*I**3 - J*J
    bad = set(factor(2*D)) | {2} if D else {2}
    for p in sorted(bad):
        r = qp_soluble(g, p)
        if r is False: return False
        if r is None: return None
    return True

def has_rational_point(g, B=300):
    for v in range(0, B+1):
        for u in (range(1, B+1) if v else range(1, 2)):
            if v and gcd(u, v) != 1: continue
            val = q_eval(g, u, v)
            if val < 0: continue
            r = isqrt(val)
            if r*r == val: return (u, v, r)
    return None
