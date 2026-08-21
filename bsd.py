"""
bsd.py -- combined rank machine.

Two independent engines, squeezed together:

  UPPER  descent.py  2-isogeny descent   ->  rank <= log2|Sel| + log2|Sel'| - 2
  LOWER  mw.py       point search + Neron-Tate pairing  ->  rank >= #independent

They search different objects.  descent.py looks for points on the quartic
torsors w^2 = d u^4 + a u^2 v^2 + (b/d) v^4; mw.py looks for points on E
itself and proves independence via the height pairing.  Neither subsumes
the other, so the combination pins ranks that neither engine pins alone.

When lower == upper the rank is EXACT.
When lower <  upper the difference is dim Sha[2] -- provided neither engine
was merely out of puff, which the report distinguishes.
"""
import descent, mw
from fractions import Fraction as F

def upper(a, b, B=600, nb=2_000_000, md=20):
    descent.NB[0], descent.MD[0] = nb, md
    r = descent.rank_bounds(a, b, B=B)
    return r['lo'], r['hi'], bool(r['undet_E'] or r['undet_F'])

def lower(a2, a4, a6, E=6, M=200, prec=7):
    c = (a2, a4, a6)
    pts = mw.search_points(c, E, M)
    ind, reg = mw.independent_subset(pts, c, n=prec)
    return len(ind), reg, ind

def analyse(a, b, label, Bd=600, E=6, M=200, prec=7, quiet=False):
    """E : y^2 = x^3 + a x^2 + b x   (needs 2-torsion for the upper bound)."""
    d_lo, d_hi, undet = upper(a, b, B=Bd)
    p_lo, reg, ind = lower(a, b, 0, E, M, prec)
    lo, hi = max(d_lo, p_lo), d_hi
    if lo > hi:                       # engines disagree -> something is wrong
        return dict(label=label, lo=lo, hi=hi, status="CONTRADICTION",
                    d=(d_lo, d_hi), p=p_lo, reg=reg, undet=undet)
    if lo == hi:      status = "EXACT"
    elif undet:       status = "local test inconclusive"
    else:             status = f"gap {hi-lo}"
    return dict(label=label, lo=lo, hi=hi, status=status,
                d=(d_lo, d_hi), p=p_lo, reg=reg, undet=undet)

def sweep(ns, Bd=600, E=6, M=200, prec=7):
    print(f"\n{'='*88}")
    print("  COMBINED SWEEP -- congruent number curves  y^2 = x^3 - n^2 x")
    print(f"{'='*88}")
    print("   n  | descent  | points | COMBINED | regulator   | status")
    print("  " + "-"*84)
    rows = []
    for n in ns:
        r = analyse(0, -n*n, f"n={n}", Bd, E, M, prec)
        rows.append((n, r))
        gain = ""
        if r['p'] > r['d'][0]: gain = " <- points beat descent"
        if r['status'] == "EXACT" and r['d'][0] < r['d'][1]:
            gain = " <- COMBINATION pinned it"
        print(f"  {n:>3} | [{r['d'][0]},{r['d'][1]}]{'*' if r['undet'] else ' '}    "
              f"| >= {r['p']}   | [{r['lo']},{r['hi']}]    "
              f"| {r['reg']:>10.5f}  | {r['status']}{gain}")
    return rows

if __name__ == "__main__":
    rows = sweep(range(1, 61))
    ex   = [n for n, r in rows if r['status'] == "EXACT"]
    gap  = [n for n, r in rows if r['status'].startswith("gap")]
    und  = [n for n, r in rows if r['undet'] and r['status'] != "EXACT"]
    bad  = [n for n, r in rows if r['status'] == "CONTRADICTION"]
    print("\n  " + "-"*84)
    print(f"  exact ranks : {len(ex)}/{len(rows)}")
    print(f"  Sha[2] gaps : {gap}")
    print(f"  inconclusive: {und}")
    print(f"  CONTRADICTIONS (bug indicator): {bad}")
