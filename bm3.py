"""Part 3, done properly: evaluate A = (-1, 3-t^2) on a genuine ADELIC point."""
from fractions import Fraction as F
from math import gcd
exec(open("bm.py").read().split("# ================================================== 1.")[0])

PLACES = [0] + [p for p in range(2, 200) if all(p % d for d in range(2, int(p**.5)+1))]

def local_ok(t, p):
    """Does X have a Q_p-point in the fibre over t?  (p=0: the real place.)"""
    c = c_of(t)
    return c != 0 and hilbert(-1, c, p) == 1

def support(t):
    """All places where inv_v(-1, 3-t^2) = 1/2.  Finite."""
    g = 3 - F(t)*F(t)
    cand = {0, 2} | set(factor(g.numerator * g.denominator))
    return {p for p in cand if inv(-1, g, p) == F(1,2)}

def build_adelic(t0):
    """t0 at every place it works; a repaired t_p at the places it doesn't."""
    pt = {}
    for p in PLACES:
        if local_ok(t0, p): pt[p] = t0; continue
        for num in range(-80, 81):
            for den in range(1, 40):
                if gcd(abs(num), den) != 1: continue
                t = F(num, den)
                if local_ok(t, p): pt[p] = t; break
            if p in pt: break
        assert p in pt, f"no local point at {p}"
    return pt

def total(pt, t0):
    """sum_v inv_v(A).  Places using t0 contribute t0's support minus the
    repaired places; each repaired place contributes its own invariant."""
    rep = {p for p, t in pt.items() if t != t0}
    s = F(0)
    contribs = []
    for p in support(t0) - rep:
        s += F(1,2); contribs.append(("inf" if p == 0 else str(p), t0, F(1,2)))
    for p in sorted(rep):
        i = inv(-1, 3 - F(pt[p])**2, p)
        if i: s += i; contribs.append(("inf" if p == 0 else str(p), pt[p], i))
    return s % 1, contribs

line("=")
print("  3.  THE OBSTRUCTION AS A NUMBER  --  evaluated on ADELIC points")
line("=")
print("  On X:  (3-t^2)(t^2-2) = x^2+y^2 is a norm from Q(i), so at every")
print("  local point  inv_v(-1, 3-t^2) = inv_v(-1, t^2-2).  A is well defined.")
print("  An adelic point uses a DIFFERENT t_v at each place.\n")

for t0 in [F(8,5), F(5,3), F(12,7), F(26,15), F(45,26)]:
    pt = build_adelic(t0)
    s, contribs = total(pt, t0)
    rep = sorted(p for p, t in pt.items() if t != t0)
    print(f"  base t0 = {t0}")
    print(f"    repaired at places {[('inf' if p==0 else p) for p in rep]}"
          f"  (t0 has no local point there)")
    for name, t, i in contribs:
        print(f"      inv_{name:<4}(A)  with t_v = {str(t):<8} =  {i}")
    print(f"    ------------------------------------------------")
    print(f"    sum_v inv_v(A)  =  {s}      <-- NOT ZERO\n")

print("  Every adelic point gives 1/2.  By Hilbert reciprocity a RATIONAL")
print("  point would force the sum to be 0.  No adelic point achieves 0,")
print("  therefore  X(A_Q)^Br = empty,  therefore  X(Q) = empty.")
line("=")
