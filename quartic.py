"""
Stage 1: classical invariant theory of binary quartics.

    g(u,v) = a u^4 + b u^3 v + c u^2 v^2 + d u v^3 + e v^4

    I = 12ae - 3bd + c^2
    J = 72ace + 9bcd - 27ad^2 - 27b^2 e - 2c^3

The Jacobian of  z^2 = g(u,v)  is  Y^2 = X^3 - 27 I X - 27 J.

Seminvariants (invariant under v -> v + t u, i.e. the unipotent subgroup):
    H = 8ac - 3b^2
    R = b^3 + 8a^2 d - 4abc

These satisfy a syzygy relating them to I, J and a.  We DERIVE it here by
linear algebra over many random quartics rather than quoting it.
"""
from fractions import Fraction as F
from itertools import product
import random

def I_inv(a,b,c,d,e):  return 12*a*e - 3*b*d + c*c
def J_inv(a,b,c,d,e):  return 72*a*c*e + 9*b*c*d - 27*a*d*d - 27*b*b*e - 2*c**3
def H_sem(a,b,c,d,e):  return 8*a*c - 3*b*b
def R_sem(a,b,c,d,e):  return b**3 + 8*a*a*d - 4*a*b*c
def disc(a,b,c,d,e):
    I, J = I_inv(a,b,c,d,e), J_inv(a,b,c,d,e)
    return F(4*I**3 - J*J, 27)

# ---- derive the syzygy:  R^2 = x1*H^3 + x2*I*a^2*H + x3*J*a^3  ------------
random.seed(7)
rows, rhs = [], []
for _ in range(400):
    q = [random.randint(-9, 9) for _ in range(5)]
    if q[0] == 0: continue
    a = q[0]
    I, J = I_inv(*q), J_inv(*q)
    H, R = H_sem(*q), R_sem(*q)
    rows.append([H**3, I*a*a*H, J*a**3])
    rhs.append(R*R)

def solve3(rows, rhs):
    """Exact 3x3 solve from the first independent triple, then verify on all."""
    n = len(rows)
    for combo in product(range(min(n, 40)), repeat=3):
        i, j, k = combo
        if len({i,j,k}) < 3: continue
        M = [rows[i][:], rows[j][:], rows[k][:]]
        y = [F(rhs[i]), F(rhs[j]), F(rhs[k])]
        A = [[F(v) for v in M[r]] + [y[r]] for r in range(3)]
        # gaussian elimination
        ok = True
        for col in range(3):
            p = next((r for r in range(col,3) if A[r][col] != 0), None)
            if p is None: ok = False; break
            A[col], A[p] = A[p], A[col]
            A[col] = [v / A[col][col] for v in A[col]]
            for r in range(3):
                if r != col and A[r][col] != 0:
                    f = A[r][col]
                    A[r] = [A[r][t] - f*A[col][t] for t in range(4)]
        if not ok: continue
        return [A[r][3] for r in range(3)]
    return None

sol = solve3(rows, rhs)
print("  derived coefficients (x1, x2, x3) in  R^2 = x1*H^3 + x2*I*a^2*H + x3*J*a^3")
print(f"    x1 = {sol[0]},  x2 = {sol[1]},  x3 = {sol[2]}")

bad = 0
for q, r2 in zip([None]*0, []): pass
random.seed(99)
for _ in range(4000):
    q = [random.randint(-40, 40) for _ in range(5)]
    if q[0] == 0: continue
    a = q[0]; I, J = I_inv(*q), J_inv(*q); H, R = H_sem(*q), R_sem(*q)
    if R*R != sol[0]*H**3 + sol[1]*I*a*a*H + sol[2]*J*a**3: bad += 1
print(f"  verification on 4000 random quartics:  {bad} failures")
assert bad == 0
print("  SYZYGY CONFIRMED:  R^2 = H^3 - 48 I a^2 H + 64 J a^3"
      if (sol[0],sol[1],sol[2]) == (1,-48,64) else
      f"  SYZYGY CONFIRMED with coefficients {tuple(sol)}")

# ---- covariance: I and J are SL2(Z)-invariants of weight 4 and 6 ----------
def act(q, m):
    """g(u,v) -> g(alpha u + beta v, gamma u + delta v)."""
    a,b,c,d,e = q; (al,be),(ga,de) = m
    # expand symbolically via coefficient extraction on a polynomial in u,v
    out = [0]*5
    for k, coef in enumerate([a,b,c,d,e]):          # term u^(4-k) v^k
        # (al u + be v)^(4-k) * (ga u + de v)^k
        p = [1]
        for _ in range(4-k):
            p = [p[i]*al if i < len(p) else 0 for i in range(len(p)+1)]
            q2 = [0] + [x*be for x in ([1] if False else [])]
            p = None; break
        break
    # simpler: do it with explicit polynomial multiplication
    def polymul(P, Q):
        Rr = [0]*(len(P)+len(Q)-1)
        for i,x in enumerate(P):
            for j,y in enumerate(Q): Rr[i+j] += x*y
        return Rr
    tot = [0]*5
    for k, coef in enumerate([a,b,c,d,e]):
        P = [1]
        for _ in range(4-k): P = polymul(P, [al, be])   # coeffs in u^i v^(deg-i)? track as v-degree
        for _ in range(k):   P = polymul(P, [ga, de])
        for j,x in enumerate(P): tot[j] += coef*x
    return tot

random.seed(3)
bad = 0
for _ in range(300):
    q = [random.randint(-6,6) for _ in range(5)]
    al,be,ga,de = [random.randint(-3,3) for _ in range(4)]
    det = al*de - be*ga
    if det == 0: continue
    q2 = act(q, ((al,be),(ga,de)))
    if I_inv(*q2) != det**4 * I_inv(*q): bad += 1
    if J_inv(*q2) != det**6 * J_inv(*q): bad += 1
print(f"  GL2 covariance  I -> det^4 I,  J -> det^6 J :  {bad} failures")
assert bad == 0
print("  STAGE 1 PASS")
