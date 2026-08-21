"""
cond.py -- the conductor, without a full Tate implementation.

  * Minimise (c4,c6) by the Kraus conditions.
  * For p >= 5 the conductor exponent is exact:
        f_p = 0                     p does not divide Delta
        f_p = 1                     p | Delta, p does not divide c4   (multiplicative)
        f_p = 2                     p | Delta, p | c4                 (additive)
  * For p = 2, 3 the exponent is bounded (f_2 <= 8, f_3 <= 5), so we
    enumerate the few candidates and let the L-function decide: the true
    conductor is the one for which the analytic continuation reproduces the
    Dirichlet series in its region of convergence.  That is an INDEPENDENT
    check, not a fit.
"""
from descent import factor

def vp(n, p):
    if n == 0: return 10**9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

def bs(a):
    a1,a2,a3,a4,a6 = a
    return (a1*a1+4*a2, 2*a4+a1*a3, a3*a3+4*a6,
            a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4)

def disc(a):
    b2,b4,b6,b8 = bs(a)
    return -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6

def c4c6(a):
    b2,b4,b6,b8 = bs(a)
    return b2*b2-24*b4, -b2**3+36*b2*b4-216*b6

def minimise(c4, c6):
    """Strip u^4,u^6 while Delta = (c4^3-c6^2)/1728 stays integral."""
    changed = True
    while changed:
        changed = False
        D = (c4**3 - c6*c6)//1728
        if D == 0: break
        for p in sorted(factor(abs(D))):
            if vp(c4,p) >= 4 and vp(c6,p) >= 6 and vp(D,p) >= 12:
                nc4, nc6 = c4//p**4, c6//p**6
                nD = (nc4**3 - nc6*nc6)
                if nD % 1728 == 0:
                    c4, c6 = nc4, nc6; changed = True; break
    return c4, c6

def candidates(a):
    """Yield candidate conductors, best-guess first."""
    c4, c6 = c4c6(a)
    c4, c6 = minimise(c4, c6)
    D = (c4**3 - c6*c6)//1728
    base = 1; ps = sorted(factor(abs(D)))
    for p in ps:
        if p < 5: continue
        base *= p**(1 if vp(c4,p) == 0 else 2)
    e2 = [0] if 2 not in ps else list(range(1,9))
    e3 = [0] if 3 not in ps else list(range(1,6))
    # heuristic ordering: small exponents first
    out = []
    for x in e2:
        for y in e3:
            out.append(base * 2**x * 3**y)
    return sorted(set(out)), (c4, c6, D)

def conductor(a, M=2200, probes=(2.0, 2.3)):
    """Pick the candidate whose continuation matches the Dirichlet series.
    Probes sit just above the convergence edge (Re s = 3/2) where the value
    genuinely depends on N; far to the right every N looks the same, which is
    why a probe at s=3 has no discriminating power."""
    import lfun
    cands, (c4, c6, D) = candidates(a)
    best = None
    for N in cands:
        try: L = lfun.Lfun(a, N, M=M)
        except Exception: continue
        for eps in (1,-1):
            errs = []
            for s in probes:
                d = L.dirichlet(s, M=M)
                v = L.L(s, eps)
                errs.append(abs(v-d)/max(abs(d),1e-9))
            err = max(errs)                      # must match at EVERY probe
            if best is None or err < best[0]: best = (err, N, eps)
    if best is None: return None, None, None
    err, N, eps = best
    return N, eps, err
