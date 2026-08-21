"""
rank.py -- 2-descent rank calculator for elliptic curves over Q.

    python rank.py a1 a2 a3 a4 a6      general Weierstrass
                                       y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6
    python rank.py -s A B              short form  y^2 = x^3 + A x + B

Requires a rational point of order 2 (equivalently, the cubic must have a
rational root).  Reports rank_lower <= rank <= rank_upper; when they differ
the difference is dim Sha[2] -- unless the point search or the local test
was the weak link, which the report says explicitly.
"""
import sys
from descent import *

def to_ab(a1, a2, a3, a4, a6):
    """Weierstrass -> y^2 = x^3 + a x^2 + b x, via a rational 2-torsion point."""
    b2 = a1*a1 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3*a3 + 4*a6
    P, Q, R = b2, 8*b4, 16*b6            # Z^2 = X^3 + P X^2 + Q X + R
    roots = [r for r in squarefree_divisors(R) + [d for d in range(-2000, 2001)]
             if r*r*r + P*r*r + Q*r + R == 0] if R else [0]
    if not roots: return None
    r = roots[0]
    a, b = 3*r + P, 3*r*r + 2*P*r + Q
    u = 2                                 # strip 4th powers: x -> u^2 x
    while True:
        for u in range(2, 40):
            if a % (u*u) == 0 and b % (u**4) == 0:
                a, b = a // (u*u), b // (u**4); break
        else: break
    return a, b

def report(a, b, B=1500, label=""):
    NB[0], MD[0] = 2_000_000, 20
    r = rank_bounds(a, b, B=B)
    lo, hi = r['lo'], r['hi']
    und = bool(r['undet_E'] or r['undet_F'])
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  E  : y^2 = x^3 + {a}x^2 + {b}x")
    print(f"  E' : Y^2 = X^3 + {r['a2']}X^2 + {r['b2']}X      (2-isogenous)")
    print(f"{'='*70}")
    print(f"  alpha (E)  image {r['found_E']}")
    print(f"             Selmer {r['selmer_E']}")
    print(f"  alpha'(E') image {r['found_F']}")
    print(f"             Selmer {r['selmer_F']}")
    print(f"  {'-'*66}")
    print(f"  rank(E/Q)  in  [{lo}, {hi}]")
    if lo == hi:
        print(f"  ==> RANK = {lo}   (exact: descent closed)")
    elif und:
        print(f"  ==> gap {hi-lo}, but local test inconclusive -- NOT a Sha claim")
    else:
        print(f"  ==> gap {hi-lo}.  Every torsor is everywhere locally solvable and")
        print(f"      has no rational point of height <= {B}.")
        print(f"      If rank = {lo} then dim Sha[2] = {hi-lo}, #Sha[2] = {2**(hi-lo)}")
    return r

if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        for lbl, w in [("Congruent number 5   (rank 1)",      (0,0,0,-25,0)),
                       ("Congruent number 17  (Sha)",         (0,0,0,-289,0)),
                       ("Congruent number 42  (Sha)",         (0,0,0,-1764,0)),
                       ("Curve 37a1  (first rank 1 curve)",   (0,0,1,-1,0)),
                       ("Curve 389a1 (first rank 2 curve)",   (0,1,1,-2,0)),
                       ("Curve 5077a1(first rank 3 curve)",   (0,0,1,-7,6))]:
            ab = to_ab(*w)
            if ab is None:
                print(f"\n  {lbl}:  no rational 2-torsion -- this method does not apply")
                continue
            report(*ab, label=lbl)
    elif args[0] == "-s":
        A, B_ = int(args[1]), int(args[2])
        ab = to_ab(0, 0, 0, A, B_)
        print("no rational 2-torsion" if ab is None else "") or report(*ab, label=f"y^2 = x^3 + {A}x + {B_}")
    else:
        ab = to_ab(*[int(x) for x in args[:5]])
        print("no rational 2-torsion" if ab is None else "") or report(*ab, label=" ".join(args[:5]))
