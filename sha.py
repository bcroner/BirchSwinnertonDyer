"""
Sha, made visible.

Two genus-1 curves with points in EVERY completion of Q (R and every Q_p)
and NO rational points whatsoever.  Each is a nontrivial element of
Sha(J/Q) for its Jacobian J.  This is precisely the gap between what
descent can prove and what is actually true.
"""
from math import isqrt, gcd

# ---------------------------------------------------------------- utilities
def vp(n, p):
    """p-adic valuation."""
    if n == 0: return 10**9
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v

def is_p_adic_square(t, p):
    """Is t a square in Q_p?  (t a nonzero integer.)"""
    if t == 0: return True
    v = vp(t, p)
    if v % 2: return False
    u = t // p**v
    if p == 2:
        return u % 8 == 1                      # unit square in Z_2 iff = 1 mod 8
    return pow(u % p, (p - 1)//2, p) == 1      # unit square in Z_p iff QR mod p

def icbrt(n):
    """Exact integer cube root, or None."""
    neg = n < 0; n = abs(n)
    r = round(n ** (1/3)) if n < 10**15 else int(n ** (1/3))
    for c in range(max(0, r-2), r+3):
        if c**3 == n: return -c if neg else c
    return None

def is_sq(n):
    if n < 0: return False
    r = isqrt(n); return r*r == n

# ================================================================ CURVE A
# Lind-Reichardt:  C_A :  z^2 = 2 x^4 - 17 w^4     (in P(1,2,1))
# ------------------------------------------------------------------------
def A_local(p, B=40):
    """Find a Q_p-point: primitive (x,w) with 2x^4-17w^4 a square in Q_p."""
    for x in range(-B, B+1):
        for w in range(-B, B+1):
            if (x, w) == (0, 0) or gcd(x, w) != 1: continue
            t = 2*x**4 - 17*w**4
            if t != 0 and is_p_adic_square(t, p):
                return (x, w, t)
    return None

def A_real():
    for x in range(1, 20):
        t = 2*x**4 - 17
        if t >= 0: return (x, 1, t)
    return None

def A_global(B):
    """Exhaustive search for a rational point of height <= B."""
    found = []
    for x in range(0, B+1):
        x4 = 2*x**4
        for w in range(1, B+1):
            if gcd(x, w) != 1: continue
            t = x4 - 17*w**4
            if t >= 0 and is_sq(t): found.append((x, w, isqrt(t)))
    # the w=0 "point at infinity" branch: z^2 = 2x^4 needs sqrt(2) rational
    return found

# ================================================================ CURVE B
# Selmer:  C_B :  3 x^3 + 4 y^3 + 5 z^3 = 0        (smooth plane cubic)
# ------------------------------------------------------------------------
def B_f(x, y, z):  return 3*x**3 + 4*y**3 + 5*z**3
def B_grad(x, y, z): return (9*x**2, 12*y**2, 15*z**2)

def B_local(p, N):
    """
    Hensel: if v_p(f(a)) > 2 * min_i v_p(d f/d x_i (a)) then a lifts to a
    genuine Q_p-point.  A primitive solution has a unit coordinate, so we
    may normalise that coordinate to 1 and scan the other two mod p^N.
    """
    M = p**N
    for slot in range(3):
        for u in range(M):
            for v in range(M):
                a = [u, v]; a.insert(slot, 1)
                x, y, z = a
                fv  = vp(B_f(x, y, z) % (M*M*M) or M*M*M, p)
                gv  = min(vp(g % M or M, p) for g in B_grad(x, y, z))
                if fv > 2*gv and 2*gv + 1 <= N:
                    return (x, y, z, fv, gv)
    return None

def B_global(B):
    """Exhaustive: 5z^3 = -(3x^3+4y^3).  Solve for z, test perfect cube."""
    found = []
    for x in range(-B, B+1):
        x3 = 3*x**3
        for y in range(-B, B+1):
            if (x, y) == (0, 0): continue
            n = -(x3 + 4*y**3)
            if n % 5: continue
            z = icbrt(n // 5)
            if z is not None and gcd(gcd(abs(x), abs(y)), abs(z)) == 1:
                found.append((x, y, z))
    return found

# ================================================================ REPORT
def line(c="-"): print(c*74)

print(); line("=")
print("  CURVE A  (Lind-Reichardt):   z^2 = 2 x^4 - 17 w^4")
line("=")
print("  Bad primes: 2, 17.  For p >= 5 with good reduction, Hasse-Weil gives")
print("  #C(F_p) >= p + 1 - 2*sqrt(p) > 0, and Hensel lifts.  So it suffices")
print("  to exhibit points at R and at the bad primes explicitly.\n")

x, w, t = A_real()
print(f"  R      :  (x,w)=({x},{w})  ->  z^2 = {t} >= 0                      OK")
for p in [2, 17, 3, 5, 7, 11, 13, 19, 23, 29, 31, 37, 41, 43, 47]:
    r = A_local(p)
    assert r, f"no Q_{p} point found"
    x, w, t = r
    print(f"  Q_{p:<5}:  (x,w)=({x},{w})  ->  z^2 = {t:<8} is a square in Q_{p}   OK")

print("\n  => C_A has points EVERYWHERE LOCALLY.\n")
BND = 900
pts = A_global(BND)
print(f"  Global search, all coprime (x,w) with 0 <= x,w <= {BND}"
      f"  [{(BND+1)**2:,} pairs]")
print(f"  Rational points found: {len(pts)}")
print("  The w=0 branch needs z^2 = 2x^4, i.e. sqrt(2) in Q.  Impossible.")
print("\n  => NO RATIONAL POINTS.   (Lind 1940 / Reichardt 1942, by descent.)")

# Jacobian via the classical quartic invariants of f = a x^4+b x^3+c x^2+d x+e
a, b, c, d, e = 2, 0, 0, 0, -17
I = 12*a*e - 3*b*d + c**2
J = 72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*b**2*e - 2*c**3
print(f"\n  Jacobian of C_A:  I = {I}, J = {J}")
print(f"    y^2 = x^3 - 27*I*x - 27*J  =  x^3 + {-27*I}x")
n = -27*I
u = 1
while n % 3**4 == 0: n //= 3**4; u *= 3
while n % 2**4 == 0: n //= 2**4; u *= 2
print(f"    minimal-ish twist:  y^2 = x^3 + {n}x        ({n} = 2^3 * 17)")
print("    C_A is a nontrivial element of Sha(J_A/Q)[2].")

print(); line("=")
print("  CURVE B  (Selmer 1951):      3 x^3 + 4 y^3 + 5 z^3 = 0")
line("=")
print("  Bad primes: 2, 3, 5.  Hensel criterion used:  v_p(f) > 2*min_i v_p(df/dx_i).\n")
print("  R      :  x=1, y=1  ->  5z^3 = -7, z real                        OK")
for p, N in [(2,6), (3,5), (5,4), (7,3), (11,2), (13,2), (17,2), (19,2), (23,2)]:
    r = B_local(p, N)
    assert r, f"no Q_{p} point found"
    x, y, z, fv, gv = r
    print(f"  Q_{p:<5}:  (x,y,z)=({x},{y},{z})  v_p(f)={fv} > 2*{gv}=2*v_p(grad)   OK")

print("\n  => C_B has points EVERYWHERE LOCALLY.\n")
BND = 600
pts = B_global(BND)
print(f"  Global search, all (x,y,z) with |x|,|y| <= {BND}"
      f"  [{(2*BND+1)**2:,} pairs]")
print(f"  Rational points found: {len(pts)}")
print("\n  => NO RATIONAL POINTS.   (Selmer 1951.)")
print("     C_B is a nontrivial element of Sha(J_B/Q)[3].")
print(); line("=")
