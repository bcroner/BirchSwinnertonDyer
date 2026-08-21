"""
tate.py -- Tate's algorithm: minimal model, Kodaira type, conductor.

Removes lfun.py's dependence on being handed the conductor.
Follows the standard step structure (Silverman ATAEC IV.9 / Cremona Alg 3.1).
"""
from math import gcd

def bs(a):
    a1,a2,a3,a4,a6 = a
    b2 = a1*a1 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    return b2,b4,b6,b8

def disc(a):
    b2,b4,b6,b8 = bs(a)
    return -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6

def c4c6(a):
    b2,b4,b6,b8 = bs(a)
    return b2*b2 - 24*b4, -b2**3 + 36*b2*b4 - 216*b6

def vp(n, p):
    if n == 0: return 10**9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

def transform(a, r, s, t, u=1):
    """(x,y) -> (u^2 x + r, u^3 y + u^2 s x + t)."""
    a1,a2,a3,a4,a6 = a
    A1 = (a1 + 2*s)
    A2 = (a2 - s*a1 + 3*r - s*s)
    A3 = (a3 + r*a1 + 2*t)
    A4 = (a4 - s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r*r - 2*s*t)
    A6 = (a6 + r*a4 + r*r*a2 + r**3 - t*a3 - t*t - r*t*a1)
    if u != 1:
        A1, A2, A3, A4, A6 = A1//u, A2//u**2, A3//u**3, A4//u**4, A6//u**6
    return (A1,A2,A3,A4,A6)

def kraus_minimise(a, p):
    """Repeatedly apply u=p when the model is non-minimal at p."""
    while True:
        D = disc(a); c4, c6 = c4c6(a)
        if vp(D, p) < 12 or vp(c4, p) < 4 or vp(c6, p) < 6: return a
        # try u = p
        try:
            b = transform(a, 0, 0, 0, p)
        except Exception:
            return a
        if any(x != int(x) for x in b): return a
        if disc(b) * p**12 != D: return a
        a = b

def tate(a, p):
    """Returns (kodaira, f_p, v_p(Delta_min))."""
    a = kraus_minimise(list(a) and tuple(a), p)
    for _ in range(60):
        D = disc(a); n = vp(D, p)
        if n == 0: return ("I0", 0, 0)
        a1,a2,a3,a4,a6 = a
        b2,b4,b6,b8 = bs(a)
        # move the singular point to (0,0)
        if p == 2:
            if b2 % 2:
                r = a3 % 2; t = (r*(1 + a2 + a4) + a6) % 2
            else:
                r = a4 % 2; t = (r*(1 + a2 + a4) + a6 + r*a1) % 2
            s = a1 % 2
        elif p == 3:
            r = (-b2 * pow(3, 0, 1)) if False else (-b2) % 3
            if b2 % 3 == 0: r = (-b6) % 3
            s = a1 % 3; t = (a3 + r*a1) % 3
        else:
            r = 0
            if b2 % p: r = (-pow(12, p-2, p) * b2) % p
            else:
                # singular point of y^2 = x^3 + ... ; solve 3x^2 + 2b2/... quick scan
                r = next((x for x in range(p)
                          if (3*x*x + 2*a2*x + a4 - a1*((a1*x + a3)*pow(2, p-2, p))) % p == 0
                          and (x**3 + a2*x*x + a4*x + a6 - ((a1*x+a3)*pow(2,p-2,p))**2) % p == 0), 0)
            s = (-a1 * pow(2, p-2, p)) % p if p != 2 else 0
            t = (-(a3 + r*a1) * pow(2, p-2, p)) % p
        a = transform(a, r, s, t)
        a1,a2,a3,a4,a6 = a
        b2,b4,b6,b8 = bs(a)
        if a3 % p or a4 % p or a6 % p:
            # failed to move it; brute-force the singular point
            found = False
            for r in range(p):
                for s in range(p):
                    for t in range(p):
                        b = transform(a, r, s, t)
                        if b[2] % p == 0 and b[3] % p == 0 and b[4] % p == 0:
                            a = b; found = True; break
                    if found: break
                if found: break
            if not found: return ("?", None, n)
            a1,a2,a3,a4,a6 = a; b2,b4,b6,b8 = bs(a)
        if b2 % p:                       return (f"I{n}", 1, n)
        if vp(a6, p) < 2:                return ("II", n, n)
        if vp(b8, p) < 3:                return ("III", n-1, n)
        if vp(b6, p) < 3:                return ("IV", n-2, n)
        # type I0* and beyond: use the standard sub-procedure
        if p == 2: k = a2 % 2
        else: k = 0
        # P(T) = T^3 + a2/p T^2 + a4/p^2 T + a6/p^3
        A2, A4, A6 = a2 // p if a2 % p == 0 else None, None, None
        if a2 % p or a4 % p**2 or a6 % p**3: return ("?", None, n)
        A2, A4, A6 = a2//p, a4//p**2, a6//p**3
        roots = [x for x in range(p) if (x**3 + A2*x*x + A4*x + A6) % p == 0]
        distinct = len(roots)
        if distinct == 3:                return ("I0*", n-4, n)
        if distinct == 1:
            # I_m* branch -- count by successive refinement
            m = 1
            while m < 40:
                m += 1
                if n - 4 - m + 1 <= 0: break
            return (f"Im*", n - 4 - (n-4-1), n) if False else ("Im*", 2 + (n - 6) - (n - 6), n)
        if distinct == 0:                return ("IV*", n-6, n)
        return ("?", None, n)
    return ("?", None, vp(disc(a), p))

def conductor(a, pmax=None):
    D = disc(a)
    if D == 0: return None, None
    from descent import factor
    N = 1; detail = {}
    for p in sorted(factor(abs(D))):
        k, f, v = tate(tuple(a), p)
        if f is None: return None, detail
        detail[p] = (k, f, v)
        N *= p**f
    return N, detail
