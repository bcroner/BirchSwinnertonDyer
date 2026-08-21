"""
Stage 2: the covariants of a binary quartic, and the map to its Jacobian.

A binary quartic g has two basic covariants beyond itself:
    H  = Hessian, degree 4
    G  = Jacobian of (g, H), degree 6
The claim to establish: for suitable constants alpha, beta the map
    X = alpha * H/g ,   Y = beta * G/g^(3/2)
sends  z^2 = g(u,v)  to  Y^2 = X^3 - 27 I X - 27 J.
Equivalently, the POLYNOMIAL identity
    beta^2 G^2 = alpha^3 H^3 - 27 I alpha H g^2 - 27 J g^3
holds identically in u,v.  We solve for alpha and beta.
"""
from fractions import Fraction as F
import random
from quartic import I_inv, J_inv

# homogeneous poly of degree n as [c_0..c_n]  meaning  sum c_k u^(n-k) v^k
def pmul(P, Q):
    R = [F(0)] * (len(P) + len(Q) - 1)
    for i, x in enumerate(P):
        for j, y in enumerate(Q): R[i+j] += F(x) * F(y)
    return R

def padd(P, Q):
    n = max(len(P), len(Q)); R = [F(0)]*n
    for i, x in enumerate(P): R[i] += F(x)
    for i, x in enumerate(Q): R[i] += F(x)
    return R

def pscale(P, s): return [F(s)*F(x) for x in P]

def d_u(P):
    """d/du of sum c_k u^(n-k) v^k  -> degree n-1."""
    n = len(P) - 1
    return [F(P[k]) * (n - k) for k in range(n)]          # k = 0..n-1

def d_v(P):
    n = len(P) - 1
    return [F(P[k]) * k for k in range(1, n + 1)]

def pzero(P): return all(F(x) == 0 for x in P)

def psub(P, Q): return padd(P, pscale(Q, -1))

def hessian(g):
    return psub(pmul(d_u(d_u(g)), d_v(d_v(g))), pmul(d_u(d_v(g)), d_u(d_v(g))))

def bigG(g, H):
    return psub(pmul(d_u(g), d_v(H)), pmul(d_v(g), d_u(H)))

# ---- find alpha, beta -----------------------------------------------------
random.seed(11)
def trial(g):
    I, J = I_inv(*g), J_inv(*g)
    H = hessian(g); G = bigG(g, H)
    gg = [F(x) for x in g]
    g2 = pmul(gg, gg); g3 = pmul(g2, gg)
    G2 = pmul(G, G); H3 = pmul(pmul(H, H), H)
    return I, J, H, G, gg, g2, g3, G2, H3

# For each candidate alpha, RHS = alpha^3 H^3 - 27 I alpha H g^2 - 27 J g^3.
# Ask whether RHS is a constant multiple of G^2, and record that constant.
cands = [F(n, d) for d in (1,2,3,4,6,9,12,18,27,36,54,108) for n in range(-12, 13) if n]
g = [random.randint(-7, 7) for _ in range(5)]
I, J, H, G, gg, g2, g3, G2, H3 = trial(g)
hits = []
for al in cands:
    RHS = psub(psub(pscale(H3, al**3), pscale(pmul(H, g2), 27*I*al)), pscale(g3, 27*J))
    # is RHS = lam * G2 ?
    lam = None; ok = True
    for r, s in zip(RHS, G2):
        if s == 0:
            if r != 0: ok = False; break
            continue
        cur = F(r) / F(s)
        if lam is None: lam = cur
        elif lam != cur: ok = False; break
    if ok and lam not in (None, 0): hits.append((al, lam))
print(f"  test quartic g = {g}   (I={I}, J={J})")
print(f"  alpha values making RHS proportional to G^2: {hits}")

al, lam = hits[0]
print(f"\n  taking alpha = {al},  beta^2 = {lam}")

# ---- verify the identity on many random quartics --------------------------
bad = 0
random.seed(21)
for _ in range(300):
    g = [random.randint(-9, 9) for _ in range(5)]
    if g[0] == 0: continue
    I, J, H, G, gg, g2, g3, G2, H3 = trial(g)
    RHS = psub(psub(pscale(H3, al**3), pscale(pmul(H, g2), 27*I*al)), pscale(g3, 27*J))
    if not pzero(psub(RHS, pscale(G2, lam))): bad += 1
print(f"  polynomial identity on 300 random quartics: {bad} failures")
assert bad == 0

# ---- end-to-end: a rational point on z^2=g gives a point on the Jacobian ---
print("\n  END-TO-END CHECK: quartic point -> Jacobian point")
def jac_point(g, u, v, z):
    """(u,v,z) with z^2=g(u,v)  ->  (X,Y) on Y^2 = X^3 - 27 I X - 27 J."""
    H = hessian(g); G = bigG(g, H)
    ev = lambda P: sum(F(c) * F(u)**(len(P)-1-k) * F(v)**k for k, c in enumerate(P))
    Hv, Gv = ev(H), ev(G)
    X = al * Hv / F(z)**2
    Y = F(lam).numerator and None
    return X, Gv, Hv

import math
ok = 0; tot = 0
random.seed(5)
for _ in range(4000):
    g = [random.randint(-6, 6) for _ in range(5)]
    if g[0] == 0: continue
    u, v = random.randint(-4, 4), random.randint(-4, 4)
    if (u, v) == (0, 0): continue
    val = sum(g[k] * u**(4-k) * v**k for k in range(5))
    if val <= 0: continue
    r = math.isqrt(val)
    if r * r != val: continue
    I, J = I_inv(*g), J_inv(*g)
    H = hessian(g); G = bigG(g, H)
    ev = lambda P: sum(F(c) * F(u)**(len(P)-1-k) * F(v)**k for k, c in enumerate(P))
    X = al * ev(H) / F(r)**2
    Y2 = X**3 - 27*I*X - 27*J
    Yc = F(lam) * ev(G)**2 / F(r)**6          # beta^2 G^2 / z^6
    tot += 1
    if Y2 == Yc: ok += 1
print(f"  quartic points mapped: {tot},  landing on the Jacobian: {ok}")
assert tot and ok == tot
print("  STAGE 2 PASS")
