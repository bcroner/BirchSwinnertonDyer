"""
zeros.py -- zeros of L(E,s) on the critical line.

Lambda(s) = N^(s/2)(2pi)^-s Gamma(s) L(E,s),  Lambda(s) = w Lambda(2-s).

On the critical line s = 1 + it the functional equation forces Lambda to be
real when w = +1 and purely imaginary when w = -1 (times a fixed phase), so
zeros are found by sign changes of a REAL function -- no argument principle
needed.

Everything is built on the same split-integral formula lfun.py already uses,
extended to complex s.  Gamma(s,x) = e^-x * int_0^inf (x+u)^(s-1) e^-u du is
computed by the same Simpson rule; (x+u)^(s-1) = exp((s-1) log(x+u)) is fine
for complex s because x+u is real and positive.
"""
import cmath, math
from lfun import a_n_table, discriminant, cinvs, primes_to

def gamma_inc_c(s, x, steps=400, cut=45.0):
    """Gamma(s,x) for complex s, real x > 0."""
    h = cut/steps; tot = 0
    lx = None
    for i in range(steps+1):
        u = i*h
        f = cmath.exp((s-1)*cmath.log(x+u)) * math.exp(-u)
        w = 1 if i in (0, steps) else (4 if i % 2 else 2)
        tot += w*f
    return math.exp(-x) * tot * h/3

class LZ:
    def __init__(self, w_coeffs, N, eps, M=1600):
        self.a = a_n_table(w_coeffs, M); self.M = M
        self.N = N; self.eps = eps; self.rootN = math.sqrt(N)
    def Lambda(self, s):
        tot = 0
        for n in range(1, self.M+1):
            an = self.a[n]
            if an == 0: continue
            x = 2*math.pi*n/self.rootN
            if x > 42: break
            t1 = cmath.exp((s/2)*math.log(self.N)) * (2*math.pi)**(-s) \
                 * gamma_inc_c(s, x) * cmath.exp(-s*math.log(n))
            t2 = cmath.exp(((2-s)/2)*math.log(self.N)) * (2*math.pi)**(-(2-s)) \
                 * gamma_inc_c(2-s, x) * cmath.exp(-(2-s)*math.log(n))
            tot += an*(t1 + self.eps*t2)
        return tot
    def F(self, t):
        """Real-valued function on the critical line whose zeros are L's zeros."""
        v = self.Lambda(complex(1.0, t))
        return v.real if self.eps == 1 else v.imag

def order_at_centre(L, ts=(0.4,0.2,0.1,0.05)):
    """F(t) ~ c t^r near t=0; successive ratios tend to 2^r."""
    vs = [L.F(t) for t in ts]
    return vs, [vs[i]/vs[i+1] for i in range(len(vs)-1) if vs[i+1] != 0]

def find_zeros(L, tmax=8.0, step=0.05):
    """Sign changes of F on (0, tmax], refined by bisection."""
    out = []
    prev_t, prev = step, L.F(step)
    t = step*2
    while t <= tmax:
        cur = L.F(t)
        if prev == 0 or (prev < 0) != (cur < 0):
            lo, hi = prev_t, t
            for _ in range(50):
                mid = (lo+hi)/2
                if (L.F(lo) < 0) != (L.F(mid) < 0): hi = mid
                else: lo = mid
            out.append((lo+hi)/2)
        prev_t, prev = t, cur
        t += step
    return out

def sample_points(N, k=4):
    """Sample points must sit INSIDE the gap to the first zero, and that gap
    shrinks like 2pi/log N.  Fixed points that work at N=37 are already past
    the first zero at N=9248, which makes the ratios erratic."""
    import math
    T = min(1.6, 6.0/math.log(N))
    return tuple(T/2**i for i in range(k))

def analytic_rank(L, ts=None, verbose=False):
    if ts is None: ts = sample_points(L.N)
    """Order of vanishing at the centre, read on the CRITICAL LINE.

    F(t) = Re/Im Lambda(1+it) is real and behaves as c*t^r near t=0, so
    F(t)/F(t/2) -> 2^r.  Working on the line rather than the real axis is far
    better conditioned: on the real axis the same estimate needed 10x the
    quadrature resolution to reach 7.94 for 5077a1, while here the raw ratios
    approach 8 from below and Richardson approaches it from above, bracketing
    the answer.
    """
    import math
    vs = [L.F(t) for t in ts]
    if any(v == 0 for v in vs): return None, vs, [], False
    rs = [vs[i]/vs[i+1] for i in range(len(vs)-1)]
    rich = [2*rs[i+1]-rs[i] for i in range(len(rs)-1)]
    est = rs[-1] if not rich else (rs[-1] + rich[-1])/2
    r = round(math.log(abs(est), 2)) if est > 0 else None
    # A trustworthy run has all ratios positive and increasing toward 2^r.
    # A sign flip means a zero was crossed (sample points too large); a
    # non-monotone tail means the smallest point hit the noise floor.  Large
    # conductors are where this bites -- the zeros crowd together and the
    # values at small t are tiny.
    trusted = (all(x > 0 for x in rs) and
               all(rs[i] <= rs[i+1] + 1e-9 for i in range(len(rs)-1)))
    return r, rs, rich, trusted
