"""
descent2.py -- GENERAL 2-descent.  No rational 2-torsion required.

Birch-Swinnerton-Dyer method: 2-coverings of E correspond to binary quartics
with invariants I = c4, J = 2*c6, taken up to Q-equivalence.

  UPPER  #(everywhere-locally-solvable classes) = |Sel_2|
         rank <= dim Sel_2 - dim E(Q)[2]
  LOWER  each class WITH a rational point gives a point on E, via
         X = -H(u,v)/(4z^2),  Y = (3/32) G(u,v)/z^3   (Stage 2)
         independence from the Neron-Tate height pairing (mw.py)

SELF-AUDIT: Sel_2 is a group, so the class count MUST be a power of 2.
A non-power-of-2 count is a PROOF that the enumeration box was too small;
the tool escalates the box and reports honestly if it never stabilises.
"""
from fractions import Fraction as F
from math import gcd
import mw
from quartic import I_inv, J_inv
from equiv import act, els, has_rational_point, q_eval
from enumerate_q import enumerate_quartics
from minim import canonical_Q
from covariant import hessian, bigG

ALPHA = F(-1, 4)
BETA  = F(3, 32)

def two_torsion_dim(w):
    """dim_F2 E(Q)[2], from rational roots of the 2-division cubic.
    In the model Z^2 = X^3 + b2 X^2 + 8 b4 X + 16 b6 the 2-torsion points are
    exactly the rational roots of that cubic.  A cubic has 0, 1 or 3 rational
    roots, giving E(Q)[2] = 0, Z/2, (Z/2)^2 and dim 0, 1, 2."""
    a1, a2, a3, a4, a6 = w
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
    P, Q, R = b2, 8*b4, 16*b6
    f = lambda r: r**3 + P*r*r + Q*r + R
    roots = set()
    if R == 0:
        roots.add(0)
        # remaining quadratic X^2 + P X + Q
        D = P*P - 4*Q
        if D >= 0:
            import math
            sq = math.isqrt(D)
            if sq*sq == D:
                for r in ((-P+sq), (-P-sq)):
                    if r % 2 == 0: roots.add(r//2)
    else:
        from descent import factor
        divs = {1}
        for p, e in factor(abs(R)).items():
            divs = {d * p**k for d in divs for k in range(e+1)}
        for d in divs:
            for r in (d, -d):
                if f(r) == 0: roots.add(r)
        # if one root found, factor it out and solve the quadratic exactly
        if len(roots) == 1:
            r0 = next(iter(roots))
            A1 = P + r0; B1 = Q + r0*A1          # X^2 + A1 X + B1
            D = A1*A1 - 4*B1
            if D >= 0:
                import math
                sq = math.isqrt(D)
                if sq*sq == D:
                    for r in ((-A1+sq), (-A1-sq)):
                        if r % 2 == 0 and f(r//2) == 0: roots.add(r//2)
    n = len(roots)
    return {0: 0, 1: 1, 3: 2}.get(n, 1 if n == 2 else 0)

def invariants(a1, a2, a3, a4, a6):
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    c4 = b2*b2 - 24*b4
    c6 = -b2**3 + 36*b2*b4 - 216*b6
    return c4, c6

def quartic_to_E(g, u, v, z):
    """Map a rational point on z^2=g(u,v) to Y^2 = X^3 - 27 I X - 27 J."""
    H = hessian([F(x) for x in g]); G = bigG([F(x) for x in g], H)
    ev = lambda P: sum(F(c)*F(u)**(len(P)-1-k)*F(v)**k for k, c in enumerate(P))
    if z == 0: return None
    X = ALPHA * ev(H) / F(z)**2
    Y = BETA  * ev(G) / F(z)**3
    return (X, Y)

def selmer(I, J, boxes=((24,80),(40,150),(60,250),(90,400))):
    """Escalate the search box until the ELS count is a power of two."""
    trail = []
    for amax, cmax in boxes:
        qs = enumerate_quartics(I, J, amax=amax, cmax=cmax)
        cls = {}
        for g in qs: cls.setdefault(canonical_Q(g), []).append(g)
        good, unk = [], []
        for k in cls:
            e = els(k, I, J)
            if e is True: good.append(k)
            elif e is None: unk.append(k)
        n = len(good)
        trail.append((amax, cmax, len(qs), len(cls), n, len(unk)))
        if n > 0 and (n & (n-1)) == 0 and not unk:
            return good, trail, True
    return good, trail, False

def analyse(w, label, expect=None):
    c4, c6 = invariants(*w)
    I, J = c4, 2*c6
    print(f"\n{'='*72}")
    print(f"  {label}")
    print(f"  y^2 + {w[0]}xy + {w[2]}y = x^3 + {w[1]}x^2 + {w[3]}x + {w[4]}")
    print(f"  c4 = {c4}, c6 = {c6}   =>   I = {I}, J = {J}")
    print(f"{'='*72}")
    good, trail, ok = selmer(I, J)
    for amax, cmax, nq, ncl, n, nu in trail:
        flag = "power of 2" if (n and (n & (n-1)) == 0) else "NOT a power of 2 -> box too small"
        print(f"    box |a|<={amax:<3} |c|<={cmax:<4}: {nq:>4} quartics -> {ncl:>2} classes"
              f" -> {n} ELS   {flag}{'  [+%d inconclusive]'%nu if nu else ''}")
    if not ok:
        print("    ==> ENUMERATION DID NOT STABILISE -- upper bound NOT trustworthy")
        return
    dim = len(good).bit_length() - 1
    t2 = two_torsion_dim(w)
    print(f"    |Sel_2| = {len(good)},  dim Sel_2 = {dim},  dim E(Q)[2] = {t2}")

    # lower bound: harvest points from classes that have rational points
    cE = (0, -27*I, -27*J)
    pts = []
    for g in good:
        p = has_rational_point(g, B=200)
        if p:
            P = quartic_to_E(g, *p)
            if P and mw.on_curve(P, cE): pts.append(P)
    pts += mw.search_points(cE, 4, 60)
    ind, reg = mw.independent_subset(pts, cE, n=7)
    lo, hi = len(ind), dim - t2       # rank <= dim Sel_2 - dim E(Q)[2]
    print(f"    points harvested from torsors: {len(pts)},  independent: {len(ind)}")
    print(f"    regulator = {reg:.8f}")
    print(f"  {'-'*68}")
    print(f"  rank in [{lo}, {hi}]" + (f"    ==> RANK = {lo}" if lo == hi else
          f"    gap {hi-lo}  =>  dim Sha[2] = {hi-lo}"))
    if expect is not None:
        print(f"  known rank = {expect}   {'PASS' if lo == hi == expect else 'MISS'}")

if __name__ == "__main__":
    analyse((0,0,1,-1,0),  "37a1   -- rank 1, no rational 2-torsion", 1)
    analyse((0,1,1,-2,0),  "389a1  -- rank 2, no rational 2-torsion", 2)
    analyse((0,0,1,-7,6),  "5077a1 -- rank 3, no rational 2-torsion", 3)
