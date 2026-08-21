"""
onetors.py -- 2-descent for E with EXACTLY ONE rational 2-torsion point.

    E : y^2 = (x - e) * q(x) ,   q(T) = T^2 + bT + c irreducible over Q

The etale algebra is L = Q x K with K = Q(theta), q(theta) = 0.  The descent
map is  P = (x,y) |-> x - theta  in K*/K*^2 , and the Q-component is forced
by  (x-e) * N(x-theta) = y^2 , so the whole descent lives in K*/K*^2 and

    Sel_2  <=  K(S,2) ,   S = primes of bad reduction (plus 2).

K(S,2) is computed exactly by ks2.py, so there is no search box here either.
"""
from fractions import Fraction as F
from math import gcd
from quadfield import Quad, squarefree_part, _is_rat_square
from ks2 import in_KS2, same_class, compute_KS2, splitting
from descent import factor

def curve_data(e, b, c):
    """E : y^2 = (x-e)(x^2+bx+c).  Returns (K, m) with theta = (-b + m sqrt d)/2."""
    D = b*b - 4*c
    assert D != 0 and not _is_rat_square(F(D)), "q must be irreducible"
    d = squarefree_part(D)
    m2 = D // d
    from math import isqrt
    m = isqrt(m2)
    assert m*m == m2 and m*m*d == D
    return Quad(d), m

def x_minus_theta(K, m, b, x):
    """x - theta = (x + b/2) - (m/2) sqrt d , returned in the RING basis."""
    A = F(2*x + b, 2)          # coefficient of 1
    B = F(-m, 2)               # coefficient of sqrt d
    if K.half:
        # 1, w=(1+sqrt d)/2  =>  A + B sqrt d = (A - B) * 1 + (2B) * w
        u, v = A - B, 2*B
    else:
        u, v = A, B
    assert u.denominator == 1 and v.denominator == 1, (u, v)
    return (int(u), int(v))

def scale_to_integral(K, a):
    """Multiply by a rational square to clear denominators (class unchanged)."""
    return a

def bad_primes(e, b, c):
    """Primes of bad reduction: divisors of disc of the cubic, plus 2."""
    # disc of (x-e)(x^2+bx+c) = (b^2-4c) * ((e^2+b e+c))^2
    D = b*b - 4*c
    r = e*e + b*e + c
    S = set(factor(abs(D))) | {2}
    if r: S |= set(factor(abs(r)))
    return sorted(S)

def points_on(e, b, c, B=300):
    """Rational points with small height, for validating the descent map."""
    from math import isqrt
    out = []
    for den in range(1, 12):
        d2 = den*den
        for num in range(-B, B+1):
            if gcd(abs(num), den) != 1: continue
            x = F(num, d2)
            val = (x - e)*(x*x + b*x + c)
            n, dd = val.numerator, val.denominator
            if n < 0: continue
            if isqrt(n)**2 != n or isqrt(dd)**2 != dd: continue
            out.append(x)
    return out

# ---------------------------------------------------------------- descent
def theta_basis(K, m, b, a):
    """Ring-basis (u,v) -> theta-basis (d1,d2) with the element = d1 + d2*theta.
    sqrt d = (2 theta + b)/m, so denominators are cleared by multiplying by
    m^2 -- a rational SQUARE, so the class in K*/K*^2 is unchanged."""
    u, v = a
    if K.half:  A, Bc = F(u) + F(v, 2), F(v, 2)     # = A + B sqrt d
    else:       A, Bc = F(u), F(v)
    d1 = A + Bc*F(b, m)
    d2 = Bc*F(2, m)
    den = 1
    for z in (d1, d2): den = den*z.denominator // gcd(den, z.denominator)
    d1, d2 = d1*den*den, d2*den*den                 # multiply by (den)^2
    assert d1.denominator == 1 and d2.denominator == 1
    return int(d1), int(d2)

def homogeneous_space(e, b, c, d1, d2, delta0):
    """Returns (conic, quad) as coefficient tuples.
       conic  A u^2 + B uv + C v^2 + w^2 = 0
       square delta0 * ( P u^2 + Q uv + R v^2 - e w^2 ) = z^2 """
    A  = d2
    B  = 2*d1 - 2*b*d2
    C  = -b*d1 + (b*b - c)*d2
    P  = d1
    Q  = -2*c*d2
    R  = b*c*d2 - c*d1
    return (A, B, C, 1), (delta0*P, delta0*Q, delta0*R, -delta0*e)

def conic_point(A, B, C, D, lim=400):
    """Rational point on A u^2 + Buv + Cv^2 + Dw^2 = 0.
    For each (u,w) this is a quadratic in v, so the search is O(lim^2).

    The scan order matters enormously.  Scanning u from -lim upward returns
    whatever ugly point comes first, and parametrising from a large point
    produces a quartic with huge coefficients whose local tests then fail.
    Spiralling u outward from 0 finds the SMALLEST point, giving a small
    quartic.  (On y^2=x(x^2-x+8) this was the difference between the conic
    point (-373,-204,3) and (1,0,3), and between a wrong answer and a right
    one.)"""
    from math import isqrt
    def spiral(n):
        yield 0
        for k in range(1, n+1): yield k; yield -k
    for w in range(0, lim+1):
        for u in spiral(lim):
            if (u, w) == (0, 0): continue
            K0 = A*u*u + D*w*w
            if C == 0:
                if B*u == 0:
                    if K0 == 0: return _prim(u, 1, w) if u or w else None
                    continue
                if K0 % (B*u): continue
                return _prim(u, -K0//(B*u), w)
            disc = B*B*u*u - 4*C*K0
            if disc < 0: continue
            r = isqrt(disc)
            if r*r != disc: continue
            for num in (-B*u + r, -B*u - r):
                if num % (2*C): continue
                return _prim(u, num//(2*C), w)
    return None

def _prim(u, v, w):
    g = gcd(gcd(abs(u), abs(v)), abs(w)) or 1
    return (u//g, v//g, w//g)

# ------------------------------------------------- exact conic solvability
import sys
sys.path.insert(0, ".")
from bm import hilbert          # tested Hilbert symbol from the Brauer-Manin work

def diagonalise(A, B, C, D):
    """A u^2 + Buv + Cv^2 + Dw^2  ->  diagonal <a1,a2,a3> modulo squares."""
    if A != 0:
        return (A, A*(4*A*C - B*B), D)
    if C != 0:
        return (C, C*(4*A*C - B*B), D)
    return (B, -B, D)                       # B uv + D w^2 ~ <1,-1,D>

def conic_solvable(A, B, C, D):
    """Exact: a ternary form <a1,a2,a3> is isotropic over Q iff
    (-a1 a2, -a1 a3)_v = 1 at every place (Hasse-Minkowski + Hilbert)."""
    a1, a2, a3 = diagonalise(A, B, C, D)
    if a1 == 0 or a2 == 0 or a3 == 0: return True     # degenerate: isotropic
    x, y = -a1*a2, -a1*a3
    places = {0, 2}
    for z in (x, y):
        places |= set(factor(abs(z)))
    return all(hilbert(x, y, p) == 1 for p in places)

def bilin(con, X, Y):
    A, B, C, D = con
    return (2*A*X[0]*Y[0] + B*(X[0]*Y[1] + X[1]*Y[0])
            + 2*C*X[1]*Y[1] + 2*D*X[2]*Y[2])

def qform(con, X):
    A, B, C, D = con
    return A*X[0]**2 + B*X[0]*X[1] + C*X[1]**2 + D*X[2]**2

def parametrise(con, P0):
    """Lines through P0 cut the conic again at  Q(D) P0 - B(P0,D) D .
    With D = s E1 + r E2 this is quadratic in (s,r), giving u,v,w as
    binary quadratic forms."""
    basis = [(1,0,0), (0,1,0), (0,0,1)]
    Es = [E for E in basis if bilin(con, P0, E) != 0 or qform(con, E) != 0]
    E1, E2 = None, None
    for i in range(len(Es)):
        for j in range(i+1, len(Es)):
            M = [P0, Es[i], Es[j]]
            det = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                 - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                 + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
            if det != 0: E1, E2 = Es[i], Es[j]; break
        if E1: break
    assert E1 is not None, "no complement found"
    out = []
    for k in range(3):
        # coefficients of s^2, s r, r^2 in  Q(D) P0[k] - B(P0,D) D[k]
        qs2 = qform(con, E1); qr2 = qform(con, E2)
        qsr = bilin(con, E1, E2)
        b1 = bilin(con, P0, E1); b2 = bilin(con, P0, E2)
        c_s2 = qs2*P0[k] - b1*E1[k]
        c_sr = qsr*P0[k] - (b1*E2[k] + b2*E1[k])
        c_r2 = qr2*P0[k] - b2*E2[k]
        out.append((c_s2, c_sr, c_r2))
    return out

def quartic_from(quad, par):
    """Substitute the parametrised (u,v,w) into P u^2+Q uv+R v^2+S w^2."""
    P, Q, R, Sc = quad
    U, V, W = par
    def mul(f, g):
        out = [0]*5
        for i, a in enumerate(f):
            for j, b in enumerate(g): out[i+j] += a*b
        return out
    def add(f, g, k):
        return [x + k*y for x, y in zip(f, g)]
    res = [0]*5
    res = add(res, mul(U, U), P)
    res = add(res, mul(U, V), Q)
    res = add(res, mul(V, V), R)
    res = add(res, mul(W, W), Sc)
    g = 0
    for x in res: g = gcd(g, abs(x))
    if g:
        sq = 1
        for p, ee in factor(g).items(): sq *= p**(ee//2*2)
        if sq > 1: res = [x//sq for x in res]
    return tuple(res)

from equiv import p_solvable as q_p_solvable, real_solvable as q_real_solvable, \
                  has_rational_point as q_has_point
from quadfield import squarefree_part as _sqf

MAXDEPTH = [18]
BUDGET = [6_000_000]

def quartic_els(g, extra_primes=()):
    """Everywhere local solvability of z^2 = g(s,r).

    The quartic is GL2(Z)-reduced first.  ELS is invariant under that action,
    but the coefficient size is not, and an unreduced quartic coming from a
    large conic point can be big enough that the p-adic tests exhaust their
    budget and return 'inconclusive' -- which then inflates Sel_2."""
    from quartic import I_inv, J_inv
    from equiv import canonical
    try:
        g = canonical(g)
    except Exception:
        pass
    if not q_real_solvable(g): return False
    I, J = I_inv(*g), J_inv(*g)
    D = 4*I**3 - J*J
    bad = {2} | set(extra_primes)
    if D: bad |= set(factor(abs(D)))
    for p in sorted(bad):
        # A branch survives to level k only while v_p(g(u0,v0)) >= k - need,
        # and for primitive (u,v) that valuation is bounded in terms of
        # v_p(disc g).  So the depth needed is driven by the discriminant,
        # not by a fixed constant.
        vd = 0
        if D:
            t = abs(D)
            while t % p == 0: t //= p; vd += 1
        depth = min(2*vd + 8, 30)
        r = q_p_solvable(g, p, maxdepth=depth, budget=BUDGET[0])
        if r is False: return False
        if r is None: return None
    return True

def curve_c4c6(e, b, c):
    """c4, c6 of E : y^2 = (x-e)(x^2+bx+c), written as y^2 = x^3+A2x^2+A4x+A6."""
    A2 = b - e
    A4 = c - e*b
    A6 = -e*c
    b2, b4, b6 = 4*A2, 2*A4, 4*A6
    return b2*b2 - 24*b4, -b2**3 + 36*b2*b4 - 216*b6

def minimise_quartic(g, c4, c6):
    """The 2-covering must have invariants exactly (I,J) = (c4, 2c6).
    Scaling g -> mu*g sends I -> mu^2 I and J -> mu^3 J, so the spurious
    factor introduced by the conic parametrisation is recoverable:

        mu = J(g)*c4 / (2*c6*I(g))

    Dividing it out gives the minimal model.  Without this the quartics carry
    v_p(disc) in the dozens and no local test can decide them."""
    from quartic import I_inv, J_inv
    from fractions import Fraction as F
    I, J = I_inv(*g), J_inv(*g)
    if I == 0 or c6 == 0: return g, None
    mu = F(J*c4, 2*c6*I)
    if mu == 0: return g, None
    num, den = mu.numerator, mu.denominator
    ng = [F(x, 1)/mu for x in g]
    if all(x.denominator == 1 for x in ng):
        gg = tuple(int(x) for x in ng)
        if I_inv(*gg) == c4 and J_inv(*gg) == 2*c6: return gg, mu
    return g, mu

def selmer_onetors(e, b, c, verbose=False):
    """dim Sel_2 for E : y^2 = (x-e)(x^2+bx+c), q irreducible."""
    K, m = curve_data(e, b, c)
    S = bad_primes(e, b, c)
    reps, bound = compute_KS2(K, S)
    sel, img, unk = [], [], []
    for a in reps:
        d1, d2 = theta_basis(K, m, b, a)
        delta0 = _sqf(K.norm(a))
        con, quad = homogeneous_space(e, b, c, d1, d2, delta0)
        if not conic_solvable(*con):
            continue                                  # fails Hasse-Minkowski
        P0 = conic_point(*con)
        if P0 is None:
            unk.append(a); sel.append(a); continue     # conservative
        par = parametrise(con, P0)
        g = quartic_from(quad, par)
        if all(x == 0 for x in g): unk.append(a); sel.append(a); continue
        C4, C6 = curve_c4c6(e, b, c)
        g, mu = minimise_quartic(g, C4, C6)
        r = quartic_els(g, extra_primes=S)
        if r is None: unk.append(a); sel.append(a); continue
        if not r: continue
        sel.append(a)
        if q_has_point(g, B=250): img.append(a)
    return K, S, reps, sel, img, unk
