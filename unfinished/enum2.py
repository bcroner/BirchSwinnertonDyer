"""
Faster quartic enumeration: loop over the seminvariant H, not over c.

    H = 8ac - 3b^2,   R = b^3 + 8a^2 d - 4abc
    27 R^2 = -H^3 + 48 I a^2 H - 64 J a^3   =: phi(H)

phi(H) >= 0 is necessary, and phi(H)/27 must be a PERFECT SQUARE.  That is
an extremely strong condition, and it is cheap to pre-filter:
  * phi(H) = 0 mod 27      kills 26/27 immediately
  * quadratic-residue sieve mod 64, 63, 65 before any isqrt

Search region: |H| <= M * a * sqrt(I).  phi's stationary points are at
H = +-4a sqrt(I), so this is the natural scale; M is the safety factor.
"""
from math import isqrt, gcd, sqrt
from quartic import I_inv, J_inv

SQ64 = {(i*i) % 64 for i in range(64)}
SQ63 = {(i*i) % 63 for i in range(63)}
SQ65 = {(i*i) % 65 for i in range(65)}

def enumerate_fast(I, J, amax, M=4.0, hard_cap=None):
    out = []
    rtI = sqrt(abs(I)) if I else 1.0
    for a in range(1, amax + 1):
        Hb = int(M * a * rtI) + 8
        if hard_cap: Hb = min(Hb, hard_cap)
        A2 = 48*I*a*a
        A3 = 64*J*a**3
        for H in range(-Hb, Hb + 1):
            phi = -H**3 + A2*H - A3
            if phi < 0 or phi % 27: continue
            q = phi // 27
            if q % 64 not in SQ64: continue
            if q % 63 not in SQ63: continue
            if q % 65 not in SQ65: continue
            R = isqrt(q)
            if R*R != q: continue
            for Rs in ({R, -R}):
                for b in range(-2*a, 2*a + 1):
                    t = H + 3*b*b
                    if t % (8*a): continue
                    c = t // (8*a)
                    n2 = Rs - b**3 + 4*a*b*c
                    if n2 % (8*a*a): continue
                    d = n2 // (8*a*a)
                    n3 = I + 3*b*d - c*c
                    if n3 % (12*a): continue
                    e = n3 // (12*a)
                    g = (a, b, c, d, e)
                    if I_inv(*g) == I and J_inv(*g) == J:
                        out.append(g)
    return sorted(set(out))
