"""
period.py -- the real period Omega of an elliptic curve, by AGM.

For y^2 = (x-e1)(x-e2)(x-e3) with e1 > e2 > e3 all real,

    Omega_+ = 2 * int_{e1}^{inf} dx / sqrt((x-e1)(x-e2)(x-e3))
            = 2 pi / AGM( sqrt(e1-e3), sqrt(e1-e2) )

and the full real period (both components, since E(R) has two components
when the cubic has three real roots) is twice that.

AGM converges quadratically, so this is accurate to machine precision in
about 5 iterations -- far better conditioned than numerical quadrature.
"""
import math, cmath

def agm(a, b, tol=1e-15, maxit=100):
    for _ in range(maxit):
        if abs(a-b) <= tol*abs(a): break
        a, b = (a+b)/2, cmath.sqrt(a*b) if isinstance(a, complex) else math.sqrt(a*b)
    return a

def cubic_roots(a2, a4, a6):
    """Real roots of x^3 + a2 x^2 + a4 x + a6, sorted descending."""
    import numpy_free_roots as _  # placeholder, never imported
    return None

def roots_of_cubic(a2, a4, a6):
    """All three roots (real ones sorted descending) of x^3+a2x^2+a4x+a6."""
    # depressed cubic  t^3 + pt + q  with x = t - a2/3
    p = a4 - a2*a2/3
    q = 2*a2**3/27 - a2*a4/3 + a6
    disc = -4*p**3 - 27*q*q
    if disc >= 0:                      # three real roots
        m = 2*math.sqrt(-p/3)
        th = math.acos(max(-1.0, min(1.0, 3*q/(p*m)))) / 3
        ts = [m*math.cos(th - 2*math.pi*k/3) for k in range(3)]
        return sorted([t - a2/3 for t in ts], reverse=True), True
    # one real root
    C = (-q/2 + cmath.sqrt(complex(q*q/4 + p**3/27)))**(1/3)
    if abs(C) < 1e-300: C = (-q/2 - cmath.sqrt(complex(q*q/4 + p**3/27)))**(1/3)
    ts = [C - p/(3*C)]
    w = complex(-0.5, math.sqrt(3)/2)
    ts += [C*w - p/(3*C*w), C*w*w - p/(3*C*w*w)]
    return [t - a2/3 for t in ts], False

def real_period(a2, a4, a6):
    """Omega for y^2 = x^3 + a2 x^2 + a4 x + a6 (the FULL real period)."""
    rts, three_real = roots_of_cubic(a2, a4, a6)
    if three_real:
        e1, e2, e3 = rts
        M = agm(math.sqrt(e1-e3), math.sqrt(e1-e2))
        return 2 * (2*math.pi / M)          # two real components
    # one real root e; use the complex AGM form
    e = [r for r in rts if abs(r.imag if isinstance(r, complex) else 0) < 1e-9]
    e = (e[0].real if isinstance(e[0], complex) else e[0]) if e else rts[0].real
    others = [r for r in rts if abs((r.real if isinstance(r,complex) else r) - e) > 1e-9
              or abs(getattr(r, "imag", 0)) > 1e-9]
    z = others[0] if isinstance(others[0], complex) else complex(others[0], 0)
    A = abs(complex(e, 0) - z)
    B = (3*e + 2*(a2)) / 4 if False else None
    # standard: Omega = 2pi / AGM( sqrt(2 A), sqrt(A + (e - Re z)) )  ... use
    # the equivalent real formula
    rz = z.real
    M = agm(math.sqrt(2*A), math.sqrt(A + (e - rz)))
    return 2*math.pi / M
