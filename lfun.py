"""
lfun.py -- the L-function of an elliptic curve, from scratch.

Built ONLY from counting points mod p:

    a_p = p + 1 - #E(F_p)                 |a_p| <= 2 sqrt(p)   (Hasse)

    L(E,s) = prod_good (1 - a_p p^-s + p^(1-2s))^-1 * prod_bad (1 - a_p p^-s)^-1
           = sum_n a_n n^-s

The product converges only for Re(s) > 3/2.  s = 1 is OUTSIDE it.  That is
the entire drama: BSD asks about a point where the definition does not reach.

Modularity (Wiles, Taylor-Wiles, BCDT) says L(E,s) is the L-function of a
weight-2 newform, hence continues to an entire function with

    Lambda(s) = N^(s/2) (2 pi)^-s Gamma(s) L(E,s) = w Lambda(2-s),  w = +-1

which gives a rapidly convergent formula valid AT s=1, via the split
    Lambda(s) = sum_n a_n/n^s * A(s,n) + w * sum_n a_n/n^(2-s) * A(2-s,n)
    A(s,n) = N^(s/2) (2 pi)^-s Gamma(s, 2 pi n / sqrt N)
"""
from math import exp, log, pi, sqrt, isqrt

# ------------------------------------------------------------- curve data
def binvs(a1,a2,a3,a4,a6):
    b2 = a1*a1+4*a2; b4 = 2*a4+a1*a3; b6 = a3*a3+4*a6
    b8 = a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    return b2,b4,b6,b8

def discriminant(w):
    b2,b4,b6,b8 = binvs(*w)
    return -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6

def cinvs(w):
    b2,b4,b6,b8 = binvs(*w)
    return b2*b2-24*b4, -b2**3+36*b2*b4-216*b6

def primes_to(n):
    sieve = bytearray([1])*(n+1); sieve[0:2] = b"\0\0"
    for i in range(2, isqrt(n)+1):
        if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(n+1) if sieve[i]]

def legendre(a,p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a,(p-1)//2,p) == 1 else -1

def a_p(w, p):
    """a_p by literally counting points on E over F_p."""
    b2,b4,b6,b8 = binvs(*w)
    D = discriminant(w); c4, c6 = cinvs(w)
    if D % p == 0:
        if c4 % p == 0: return 0                    # additive reduction
        return legendre(-c6, p)                     # split(+1)/nonsplit(-1)
    if p == 2:
        cnt = 0
        for x in range(2):
            for y in range(2):
                if (y*y + w[0]*x*y + w[2]*y - (x**3 + w[1]*x*x + w[3]*x + w[4])) % 2 == 0:
                    cnt += 1
        return 2 + 1 - (cnt + 1)
    s = 0
    for x in range(p):
        s += legendre((4*x**3 + b2*x*x + 2*b4*x + b6) % p, p)
    return -s

def a_n_table(w, M):
    """Multiplicative extension of a_p to all n <= M."""
    D = discriminant(w)
    ps = primes_to(M)
    a = [0]*(M+1); a[1] = 1
    ap = {p: a_p(w,p) for p in ps}
    for p in ps:
        pk = p; prev2, prev1 = 1, ap[p]           # a_1, a_p
        pw = [1, ap[p]]
        while pk*p <= M:
            pk *= p
            nxt = ap[p]*pw[-1] - (0 if D % p == 0 else p*pw[-2])
            pw.append(nxt)
        # fill
        for e in range(1, len(pw)):
            q = p**e
            if q > M: break
            a[q] = pw[e]
    for n in range(2, M+1):
        if a[n]: continue
        for p in ps:
            if p*p > n: break
            if n % p == 0:
                q = p
                while n % (q*p) == 0: q *= p
                a[n] = a[q]*a[n//q]
                break
        else:
            pass
    # any remaining n are prime
    for n in range(2, M+1):
        if a[n] == 0 and n in set(ps): a[n] = ap[n]
    return a

# --------------------------------------------------- incomplete gamma
def gamma_inc(s, x, steps=400, cut=45.0):
    """Gamma(s,x) = int_x^inf t^(s-1) e^-t dt = e^-x int_0^inf (x+u)^(s-1) e^-u du
    by Simpson on u in [0, cut]."""
    h = cut/steps; tot = 0.0
    for i in range(steps+1):
        u = i*h
        f = (x+u)**(s-1) * exp(-u)
        wgt = 1 if i in (0, steps) else (4 if i % 2 else 2)
        tot += wgt*f
    return exp(-x) * tot * h/3

# ----------------------------------------------------------- the L-function
class Lfun:
    def __init__(self, w_coeffs, N, M=2500):
        self.w = w_coeffs; self.N = N
        self.a = a_n_table(w_coeffs, M); self.M = M
        self.rootN = sqrt(N)
    def Lambda(self, s, eps):
        tot = 0.0
        for n in range(1, self.M+1):
            an = self.a[n]
            if an == 0: continue
            x = 2*pi*n/self.rootN
            if x > 42: break
            t1 = self.N**(s/2) * (2*pi)**(-s) * gamma_inc(s, x) / n**s
            t2 = self.N**((2-s)/2) * (2*pi)**(-(2-s)) * gamma_inc(2-s, x) / n**(2-s)
            tot += an*(t1 + eps*t2)
        return tot
    def L(self, s, eps):
        from math import gamma
        return self.Lambda(s, eps) * (2*pi)**s * self.N**(-s/2) / gamma(s)
    def dirichlet(self, s, M=None):
        M = M or self.M
        return sum(self.a[n]/n**s for n in range(1, M+1))
