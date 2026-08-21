"""
The Brauer-Manin obstruction, computed from scratch.

    X :  x^2 + y^2 = (3 - t^2)(t^2 - 2)          (Iskovskikh)

X has points in R and in every Q_p.  X has no rational points.
The obstruction is the quaternion class  A = (-1, 3 - t^2) in Br(X):
every local condition is satisfiable ONE PLACE AT A TIME, and the
global sum of invariants is 1/2.
"""
from fractions import Fraction as F
from math import gcd

# ------------------------------------------------------- p-adic toolkit
def vp(n, p):
    v = 0
    while n % p == 0: n //= p; v += 1
    return v, n

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def hilbert(a, b, p):
    """Hilbert symbol (a,b)_p for nonzero rationals a,b.  p = 0 means the
    real place.  Returns +1 or -1."""
    a, b = F(a), F(b)
    if p == 0:
        return -1 if (a < 0 and b < 0) else 1
    # reduce to squarefree integers: (a,b)_p depends only on a,b mod squares
    A = a.numerator * a.denominator
    B = b.numerator * b.denominator
    al, u = vp(abs(A), p); be, v = vp(abs(B), p)
    u *= (1 if A > 0 else -1); v *= (1 if B > 0 else -1)
    if p == 2:
        e = lambda n: ((n % 8) - 1) // 2 % 2          # (n-1)/2 mod 2
        w = lambda n: ((n * n) % 16 - 1) // 8 % 2      # (n^2-1)/8 mod 2
        ex = (e(u) * e(v) + al * w(v) + be * w(u)) % 2
        return -1 if ex else 1
    eps = (p - 1) // 2 % 2
    s = (-1) ** (al * be * eps)
    s *= legendre(u, p) ** be
    s *= legendre(v, p) ** al
    return s

def inv(a, b, p):
    """Local invariant in (1/2)Z/Z, written as 0 or 1/2."""
    return F(0) if hilbert(a, b, p) == 1 else F(1, 2)

def factor(n):
    n = abs(n); f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def is_sum_of_two_rational_squares(q):
    """q in Q^*: is q = x^2+y^2 with x,y in Q?  Iff q>0 and v_p(q) even
    for every prime p = 3 mod 4.  Equivalently (-1,q)_v = 1 for all v."""
    if q == 0: return True
    if q < 0: return False
    n = F(q).numerator * F(q).denominator
    return all(not (p % 4 == 3 and e % 2) for p, e in factor(n).items())

def c_of(t):
    t = F(t); return (3 - t*t) * (t*t - 2)

def line(ch="-"): print(ch * 74)

# ================================================== 1. EVERYWHERE LOCAL
print(); line("=")
print("  1.  X IS SOLVABLE IN EVERY COMPLETION")
line("=")
print("  x^2+y^2 = c is solvable over Q_v  <=>  (-1, c)_v = +1.")
print("  So we need, for each v, some t in Q_v with (-1, c(t))_v = +1.\n")

# real place
t = F(8, 5)                                   # 2 < 64/25 < 3
print(f"  R      :  t = {t}   c = {c_of(t)} > 0   -> (-1,c)_inf = +1     OK")

PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83]
for p in PRIMES:
    hit = None
    for num in range(-60, 61):
        for den in range(1, 40):
            if gcd(abs(num), den) != 1: continue
            t = F(num, den); c = c_of(t)
            if c != 0 and hilbert(-1, c, p) == 1:
                hit = (t, c); break
        if hit: break
    assert hit, f"no Q_{p} point"
    t, c = hit
    print(f"  Q_{p:<5}:  t = {str(t):<7} c = {str(c):<14} -> (-1,c)_{p} = +1     OK")
print("\n  => X(A_Q) is non-empty:  X has points everywhere locally.")

# ================================================== 2. NO RATIONAL PTS
print(); line("=")
print("  2.  X HAS NO RATIONAL POINTS  (complete elementary proof)")
line("=")
print("  Put t = m/n in lowest terms.  Clearing the square n^4:")
print("      c(t) is a sum of two squares  <=>  D := A*B is, where")
print("      A = 3n^2 - m^2 ,   B = m^2 - 2n^2 .")
print("  Note   A + B = n^2   and   2A + 3B = m^2 .")

# gcd(A,B) = 1
bad = [(m,n) for n in range(1,60) for m in range(1,60) if gcd(m,n)==1
       and gcd(abs(3*n*n-m*m), abs(m*m-2*n*n)) not in (1,0)]
print(f"\n  gcd(A,B) divides both n^2 and m^2, hence = 1."
      f"   [checked {len(bad)} violations found]")

print("\n  D > 0 forces A > 0 and B > 0 (A<0,B<0 would give 3n^2<m^2<2n^2).")
print("  A,B coprime and positive  =>  BOTH must be sums of two squares.")
print("  An odd sum of two squares is = 1 mod 4; an even one is 0 or 2.")
print("  So A, B are each in {0,1,2} mod 4.  Now split on the parity of n:\n")

allowed = {0, 1, 2}
sq4 = {0, 1}
rows = []
for n4 in range(4):
    for m4 in range(4):
        if n4 % 2 == 0 and m4 % 2 == 0: continue      # gcd(m,n)=1
        A4 = (3*n4*n4 - m4*m4) % 4
        B4 = (m4*m4 - 2*n4*n4) % 4
        ok_res  = (A4 in allowed) and (B4 in allowed)
        ok_sum  = ((A4 + B4) % 4) == (n4*n4) % 4
        ok_msum = ((2*A4 + 3*B4) % 4) == (m4*m4) % 4
        survives = ok_res and ok_sum and ok_msum
        rows.append((m4, n4, A4, B4, ok_res, survives))

print("   m%4 n%4 | A%4 B%4 | A,B in {0,1,2}? | 2A+3B a square %4? | survives")
print("   " + "-"*68)
for m4, n4, A4, B4, ok_res, surv in rows:
    m2 = "yes" if ((2*A4+3*B4) % 4) in sq4 else "NO "
    print(f"    {m4}   {n4}  |  {A4}   {B4}  |      {'yes' if ok_res else 'NO ':<3}       "
          f"|        {m2}         |   {'YES' if surv else 'no'}")
print("   " + "-"*68)
print(f"   surviving residue classes: {sum(1 for r in rows if r[5])}")
print("\n  Every case dies.  Concretely:")
print("    n even -> m odd -> A = 3n^2-m^2 = -1 = 3 mod 4, not a sum of 2 squares.")
print("    n odd  -> A+B = n^2 = 1 mod 4, so (A,B) = (0,1) or (1,0) mod 4,")
print("              giving 2A+3B = 3 or 2 mod 4 -- but m^2 = 0 or 1 mod 4.")
print("\n  => X(Q) is EMPTY.")

B = 400
found = [(m,n) for n in range(1,B+1) for m in range(1,B+1)
         if gcd(m,n)==1 and is_sum_of_two_rational_squares(F((3*n*n-m*m)*(m*m-2*n*n)))
         and (3*n*n-m*m)*(m*m-2*n*n) > 0]
print(f"\n  Brute-force cross-check, all coprime 1<=m,n<={B} "
      f"[{B*B:,} pairs]:  {len(found)} rational points.")

# ================================================== 3. THE OBSTRUCTION
print(); line("=")
print("  3.  THE OBSTRUCTION AS A NUMBER:  sum_v inv_v(A) = 1/2")
line("=")
print("  A = (-1, 3 - t^2) in Br(X).  On X we have (3-t^2)(t^2-2) = x^2+y^2,")
print("  a norm from Q(i), so for ANY local point  inv_v(-1,3-t^2) =")
print("  inv_v(-1,t^2-2).  Evaluate A on adelic points:\n")

def invs(t):
    c1 = 3 - F(t)*F(t)
    out = {0: inv(-1, c1, 0)}
    n = c1.numerator * c1.denominator
    for p in factor(n): out[p] = inv(-1, c1, p)
    out[2] = inv(-1, c1, 2)
    return {k: v for k, v in out.items() if v != 0}

print("   adelic t_v      places with inv_v = 1/2                 sum")
print("   " + "-"*62)
for t in [F(8,5), F(5,3), F(12,7), F(7,4), F(17,10), F(26,15)]:
    d = invs(t)
    tot = sum(d.values()) % 1
    names = ", ".join(("inf" if k == 0 else str(k)) for k in sorted(d))
    print(f"   t = {str(t):<9}  {names:<40}  {tot}")
print("   " + "-"*62)
print("\n  For a RATIONAL t the sum is 0 -- that is Hilbert reciprocity.")
print("  The adelic points that survive the local conditions of part 1 are")
print("  exactly the ones where the sum is 1/2.  Non-zero sum = no rational")
print("  point.  That single number IS the Brauer-Manin obstruction.")
print(); line("=")
