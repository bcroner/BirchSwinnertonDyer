"""
Stage 3: enumerate integer binary quartics with prescribed invariants (I,J),
up to GL2(Z)-equivalence.

Uses the Stage-1 syzygy to avoid a 4-deep loop:

    H = 8ac - 3b^2 ,  R = b^3 + 8a^2 d - 4abc
    27 R^2 = -H^3 + 48 I a^2 H - 64 J a^3

so looping over (a,b,c) fixes H, hence R^2, hence R, hence d; then
    e = (I + 3bd - c^2) / (12a)
and J is re-checked as an independent consistency test.
"""
from math import isqrt, gcd
from quartic import I_inv, J_inv, H_sem, R_sem, disc

def enumerate_quartics(I, J, amax=24, cmax=80):
    """All integer quartics with these invariants and |a|<=amax, |c|<=cmax."""
    out = []
    for a in range(1, amax + 1):
        for b in range(-2*a, 2*a + 1):
            for c in range(-cmax, cmax + 1):
                H = 8*a*c - 3*b*b
                rhs = -H**3 + 48*I*a*a*H - 64*J*a**3
                if rhs < 0 or rhs % 27: continue
                q = rhs // 27
                R = isqrt(q)
                if R*R != q: continue
                for Rs in ({R, -R}):
                    num = Rs - b**3 + 4*a*b*c
                    if num % (8*a*a): continue
                    d = num // (8*a*a)
                    num2 = I + 3*b*d - c*c
                    if num2 % (12*a): continue
                    e = num2 // (12*a)
                    g = (a, b, c, d, e)
                    if I_inv(*g) == I and J_inv(*g) == J:
                        out.append(g)
    return sorted(set(out))

# ------------------------------------------------------- GL2(Z) reduction
def translate(g, t):
    """u -> u + t v."""
    a, b, c, d, e = g
    return (a,
            b + 4*a*t,
            c + 3*b*t + 6*a*t*t,
            d + 2*c*t + 3*b*t*t + 4*a*t**3,
            e + d*t + c*t*t + b*t**3 + a*t**4)

def reverse(g):
    """u <-> v."""
    a, b, c, d, e = g
    return (e, d, c, b, a)

def reduce_quartic(g, rounds=60):
    """Crude but deterministic reduction: shrink |b| by translation, then
    flip if the trailing coefficient is smaller.  Gives a canonical-ish
    representative that we use only for de-duplication."""
    best = g
    for _ in range(rounds):
        a, b, c, d, e = best
        if a == 0: break
        t = -round(b / (4*a))
        cand = translate(best, t) if t else best
        if abs(cand[0]) > abs(cand[4]) or (abs(cand[0]) == abs(cand[4]) and cand[4] < cand[0]):
            cand = reverse(cand)
        key = lambda q: (abs(q[0]), abs(q[1]), abs(q[2]), abs(q[3]), abs(q[4]))
        if key(cand) < key(best): best = cand
        else: break
    return best

def orbit_key(g, span=6):
    """Canonical key: smallest reduced form over a small GL2(Z) neighbourhood."""
    cands = set()
    for t in range(-span, span + 1):
        h = translate(g, t)
        cands.add(reduce_quartic(h))
        cands.add(reduce_quartic(reverse(h)))
        cands.add(reduce_quartic(tuple(-x for x in h)))
    return min(cands, key=lambda q: (tuple(abs(x) for x in q), q))

def classes(I, J, amax=24, cmax=80):
    qs = enumerate_quartics(I, J, amax, cmax)
    seen = {}
    for g in qs:
        seen.setdefault(orbit_key(g), []).append(g)
    return qs, seen

if __name__ == "__main__":
    for lbl, (I, J), expect in [("37a1",   (48,   -432),    1),
                                ("389a1",  (112,  -1712),   2),
                                ("5077a1", (336,  -10800),  3)]:
        qs, cls = classes(I, J)
        print(f"\n  {lbl}:  I={I}, J={J}")
        print(f"    integer quartics found : {len(qs)}")
        print(f"    GL2(Z) classes         : {len(cls)}")
        print(f"    dim Sel_2 should be    : {expect}   (=> {2**expect} classes if all ELS)")
        for k, v in sorted(cls.items())[:12]:
            print(f"      class {k}   disc={disc(*k)}   ({len(v)} reps)")
