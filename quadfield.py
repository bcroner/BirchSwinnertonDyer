"""
quadfield.py -- arithmetic in K = Q(sqrt d), d squarefree.

    O_K = Z[w],   w = sqrt d              if d = 2,3 mod 4   (disc 4d)
                  w = (1+sqrt d)/2        if d = 1 mod 4     (disc d)

Elements are (u,v) meaning u + v*w.  Everything below is exact integer
arithmetic -- no floats in the algebra.
"""
from math import isqrt, gcd

def squarefree_part(n):
    if n == 0: return 0
    s, sign = 1, (1 if n > 0 else -1)
    n = abs(n); p = 2
    while p*p <= n:
        e = 0
        while n % p == 0: n //= p; e += 1
        if e % 2: s *= p
        p += 1
    return sign * s * n

def _is_rat_square(q):
    from fractions import Fraction as F
    q = F(q)
    if q < 0: return False
    n, dd = q.numerator, q.denominator
    return isqrt(n)**2 == n and isqrt(dd)**2 == dd

def _rat_sqrt(q):
    from fractions import Fraction as F
    q = F(q)
    return F(isqrt(q.numerator), isqrt(q.denominator))

class Quad:
    def __init__(self, d):
        d = squarefree_part(d)
        assert d not in (0, 1), "not a quadratic field"
        self.d = d
        self.half = (d % 4 == 1)
        self.disc = d if self.half else 4*d

    def norm(self, a):
        u, v = a
        if self.half:  return u*u + u*v + v*v*(1-self.d)//4
        return u*u - self.d*v*v

    def trace(self, a):
        u, v = a
        return 2*u + v if self.half else 2*u

    def mul(self, a, b):
        (u1,v1),(u2,v2) = a,b
        if self.half:
            # w^2 = w + (d-1)/4
            k = (self.d-1)//4
            return (u1*u2 + v1*v2*k, u1*v2 + u2*v1 + v1*v2)
        return (u1*u2 + self.d*v1*v2, u1*v2 + u2*v1)

    def conj(self, a):
        u, v = a
        return (u+v, -v) if self.half else (u, -v)

    def to_sqrt_basis(self, a):
        """(u,v) in the ring basis -> (A,B) with the element = A + B*sqrt d."""
        from fractions import Fraction as F
        u, v = a
        if self.half: return (F(u) + F(v, 2), F(v, 2))
        return (F(u), F(v))

    def is_square_in_K(self, a):
        """Is a a square in K*?  Write a = A + B sqrt d and a = (X + Y sqrt d)^2.
        Then A = X^2 + d Y^2 and N(a) = (X^2 - d Y^2)^2, so with n = +-sqrt N(a)

            X^2 = (A + n)/2 ,      Y^2 = (A - n)/(2d)

        and BOTH must be squares of rationals.  Checking only that a
        consistent (trace, norm) pair exists is NOT enough -- that accepts
        -2 in Q(i), whose square root i*sqrt2 does not lie in the field."""
        from fractions import Fraction as F
        if a == (0, 0): return True
        A, B = self.to_sqrt_basis(a)
        N = A*A - self.d*B*B
        if N < 0: return False
        if not _is_rat_square(N): return False
        n0 = _rat_sqrt(N)
        for n in (n0, -n0):
            X2 = (A + n) / 2
            Y2 = (A - n) / (2*self.d)
            if X2 < 0 or Y2 < 0: continue
            if _is_rat_square(X2) and _is_rat_square(Y2):
                # cross-check: 2XY must equal B
                X, Y = _rat_sqrt(X2), _rat_sqrt(Y2)
                if 2*X*Y == B or 2*X*(-Y) == B: return True
        return False

    def minkowski(self):
        """Minkowski bound; for imaginary quadratic (2/pi) sqrt|disc|."""
        from math import pi, sqrt
        if self.d < 0: return (2/pi)*sqrt(abs(self.disc))
        return sqrt(self.disc)/2

    def elements_of_norm_upto(self, B):
        """All (u,v) with 0 < |N| <= B.  For d<0 the form is positive
        definite, so this region is genuinely finite and the enumeration
        is complete -- no heuristic box."""
        assert self.d < 0, "indefinite norm form: infinitely many, needs units"
        out = []
        # bound v: N >= (something) v^2
        vmax = int(isqrt(int(4*B/abs(self.disc)))) + 2
        for v in range(-vmax, vmax+1):
            umax = int(isqrt(B + abs(self.d)*v*v)) + 2
            for u in range(-umax, umax+1):
                if (u,v) == (0,0): continue
                N = self.norm((u,v))
                if 0 < abs(N) <= B: out.append((u,v))
        return out
