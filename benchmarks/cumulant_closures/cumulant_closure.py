"""
cumulant_closure.py -- general moment hierarchy with cumulant truncation for bosonic
Lindblad lattices (revision track "Cumulants", K=3 Kerr ring).

Normal-ordered monomials on K sites are keys (m, n) with m, n tuples of length K:
    prod_i (a_i^dag)^{m_i}  prod_i a_i^{n_i}.
Operators are dicts {key: coefficient}.  d<O>/dt = <L^dag O> is generated for every
monomial of total order <= N; moments of order > N appearing on the right-hand side are
replaced by their cumulant-expansion expressions with all cumulants of order > N set to
zero (moment/cumulant partition formula, applied recursively until only moments of
order <= N remain).  Nothing about translation symmetry is assumed.
"""
import itertools, time
from math import comb, factorial
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ------------------------------------------------------------------ algebra ---
def order(key):
    return sum(key[0]) + sum(key[1])

def add(A, B, c=1.0):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) + c*v
    return {k: v for k, v in out.items() if v != 0}

def scale(c, A):
    return {k: c*v for k, v in A.items()} if c != 0 else {}

def mul(A, B):
    """Normal-ordered product.  Site by site  a^n a^dag^m = sum_k k! C(n,k) C(m,k) a^dag^{m-k} a^{n-k}."""
    out = {}
    for (mA, nA), x in A.items():
        for (mB, nB), y in B.items():
            K = len(mA)
            ranges = [range(min(nA[i], mB[i]) + 1) for i in range(K)]
            for ks in itertools.product(*ranges):
                w = x*y
                m = list(mA); n = list(nB)
                for i, k in enumerate(ks):
                    w *= factorial(k)*comb(nA[i], k)*comb(mB[i], k)
                    m[i] += mB[i] - k
                    n[i] += nA[i] - k
                kk = (tuple(m), tuple(n))
                out[kk] = out.get(kk, 0) + w
    return {k: v for k, v in out.items() if v != 0}

def comm(A, B):
    return add(mul(A, B), mul(B, A), -1.0)

def dagger(A):
    return {(n, m): np.conj(v) for (m, n), v in A.items()}

def unit(K, i):
    e = [0]*K; e[i] = 1
    return tuple(e)

def zero(K):
    return tuple([0]*K)

def a_op(K, i, c=1.0):
    return {(zero(K), unit(K, i)): c}

def ad_op(K, i, c=1.0):
    return {(unit(K, i), zero(K)): c}

def lindblad_adjoint(H, jumps, O):
    """<L^dag O> = i<[H,O]> + sum_j rate_j (<c^dag O c> - 1/2 <{c^dag c, O}>)."""
    out = scale(1j, comm(H, O))
    for rate, cop in jumps:
        cd = dagger(cop); cdc = mul(cd, cop)
        out = add(out, mul(mul(cd, O), cop), rate)
        out = add(out, mul(cdc, O), -rate/2)
        out = add(out, mul(O, cdc), -rate/2)
    return out

def kerr_lattice(K, kappa, Delta, J, eta, U, bonds):
    """H = sum_j [Delta_j n_j + U/2 a_j^dag^2 a_j^2 + i(eta_j a_j^dag - eta_j^* a_j)] - J sum_bonds (a_j^dag a_l + h.c.)"""
    Delta = np.asarray(Delta, float)*np.ones(K)
    eta = np.asarray(eta, complex)*np.ones(K)
    H = {}
    for j in range(K):
        a, ad = a_op(K, j), ad_op(K, j)
        H = add(H, mul(ad, a), Delta[j])
        H = add(H, mul(mul(ad, ad), mul(a, a)), U/2)
        H = add(H, ad, 1j*eta[j]); H = add(H, a, -1j*np.conj(eta[j]))
    for (j, l) in bonds:
        H = add(H, mul(ad_op(K, j), a_op(K, l)), -J)
        H = add(H, mul(ad_op(K, l), a_op(K, j)), -J)
    jumps = [(kappa, a_op(K, j)) for j in range(K)]
    return H, jumps

def ring_bonds(K, ring=True):
    if K == 1: return []
    if K == 2: return [(0, 1)]
    return [(k, (k+1) % K) for k in range(K)] if ring else [(k, k+1) for k in range(K-1)]

# --------------------------------------------------------- set partitions ---
def set_partitions(n):
    """All set partitions of range(n) as lists of blocks (lists of indices)."""
    if n == 0:
        yield []
        return
    for part in set_partitions(n-1):
        # put element n-1 in an existing block or a new one
        for i in range(len(part)):
            yield part[:i] + [part[i] + [n-1]] + part[i+1:]
        yield part + [[n-1]]

def all_keys(K, N, nmin=1):
    """All monomial keys with nmin <= order <= N, sorted by order then lexicographically."""
    out = []
    for total in range(nmin, N+1):
        for comp in itertools.product(range(total+1), repeat=2*K):
            if sum(comp) == total:
                out.append((tuple(comp[:K]), tuple(comp[K:])))
    return out

# --------------------------------------------------------------- closure -----
class CumulantClosure:
    def __init__(self, K, N, H, jumps, obs_keys=(), verbose=False):
        t0 = time.time()
        self.K, self.N = K, N
        self.keys = all_keys(K, N)
        self.Nv = len(self.keys)
        self.idx = {k: i+1 for i, k in enumerate(self.keys)}      # 0 is the constant key
        const = (zero(K), zero(K))
        self.idx[const] = 0
        # generate d<O>/dt for every variable
        rows, cols, vals = [], [], []
        high = {}
        def col_of(key):
            if key in self.idx: return self.idx[key]
            if key not in high: high[key] = len(high)
            return None
        pending = []
        for i, key in enumerate(self.keys):
            expr = lindblad_adjoint(H, jumps, {key: 1.0})
            for k2, v in expr.items():
                c = col_of(k2)
                if c is None: pending.append((i, k2, v))
                else: rows.append(i); cols.append(c); vals.append(v)
        for key in obs_keys:
            col_of(key)
        self.high_keys = sorted(high, key=high.get)
        self.Nh = len(self.high_keys)
        for (i, k2, v) in pending:
            rows.append(i); cols.append(1 + self.Nv + high[k2]); vals.append(v)
        self.A = sp.csr_matrix((np.array(vals, complex), (rows, cols)), shape=(self.Nv, 1 + self.Nv + self.Nh))
        self.t_generate = time.time() - t0
        # cumulant-expansion polynomials of the high moments in the variables
        t0 = time.time()
        self._memo = {}
        terms_c, terms_i, terms_t = [], [], []
        maxdeg = 1
        for j, key in enumerate(self.high_keys):
            poly = self._expand(key)
            for mono, coef in poly.items():
                terms_c.append(coef); terms_i.append([self.idx[k] for k in mono]); terms_t.append(j)
                maxdeg = max(maxdeg, len(mono))
        self.maxdeg = maxdeg
        self.term_coef = np.array(terms_c, complex)
        self.term_idx = np.array([row + [0]*(maxdeg - len(row)) for row in terms_i], int)
        self.term_tgt = np.array(terms_t, int)
        self.Nterms = len(terms_c)
        self.t_closure = time.time() - t0
        if verbose:
            print(f"K={K} N={N}: {self.Nv} moment variables, {self.Nh} closed moments, "
                  f"{self.Nterms} polynomial terms, generate {self.t_generate:.2f}s closure {self.t_closure:.2f}s")

    # -- symbolic expansion: moment of order p > N with cumulants > N set to zero --
    def _factors(self, key):
        m, n = key; f = []
        for i in range(self.K):
            f += [('d', i)]*m[i]
        for i in range(self.K):
            f += [('a', i)]*n[i]
        return f

    def _key_of(self, factors):
        m = [0]*self.K; n = [0]*self.K
        for t, i in factors:
            if t == 'd': m[i] += 1
            else: n[i] += 1
        return (tuple(m), tuple(n))

    def _poly(self, key):
        if order(key) <= self.N:
            return {(key,): 1.0}
        return self._expand(key)

    @staticmethod
    def _pmul(P, Q):
        out = {}
        for k1, v1 in P.items():
            for k2, v2 in Q.items():
                k = tuple(sorted(k1 + k2)); out[k] = out.get(k, 0) + v1*v2
        return out

    def _expand(self, key):
        """<X_1...X_p> = - sum_{partitions pi, |pi|>=2} (-1)^{|pi|-1} (|pi|-1)! prod_B <X_B>  (kappa_p = 0)."""
        if key in self._memo: return self._memo[key]
        f = self._factors(key); p = len(f)
        out = {}
        for part in set_partitions(p):
            b = len(part)
            if b < 2: continue
            coef = -((-1)**(b-1))*factorial(b-1)
            poly = {(): coef}
            for blk in part:
                poly = self._pmul(poly, self._poly(self._key_of([f[i] for i in blk])))
            for k, v in poly.items():
                out[k] = out.get(k, 0) + v
        out = {k: v for k, v in out.items() if abs(v) > 1e-14}
        self._memo[key] = out
        return out

    # -- numerics --
    def high_moments(self, x):
        vals = np.concatenate([[1.0+0j], x])
        prod = np.prod(vals[self.term_idx], axis=1)*self.term_coef
        h = np.zeros(self.Nh, complex)
        np.add.at(h, self.term_tgt, prod)
        return h

    def rhs_complex(self, x):
        vals = np.concatenate([[1.0+0j], x, self.high_moments(x)])
        return self.A @ vals

    def rhs_real(self, y):
        x = y[:self.Nv] + 1j*y[self.Nv:]
        d = self.rhs_complex(x)
        return np.concatenate([d.real, d.imag])

    def moment(self, x, key):
        key = (tuple(int(v) for v in key[0]), tuple(int(v) for v in key[1]))
        if key in self.idx:
            i = self.idx[key]
            return 1.0 if i == 0 else x[i-1]
        if key in self.high_keys:
            return self.high_moments(x)[self.high_keys.index(key)]
        # any other high moment: evaluate its closure polynomial on the fly
        vals = np.concatenate([[1.0+0j], x]); tot = 0j
        for mono, coef in self._expand(key).items():
            tot += coef*np.prod([vals[self.idx[k]] for k in mono])
        return tot

    def coherent_state(self, alpha):
        """Moments of a product coherent state (all cumulants beyond first vanish)."""
        alpha = np.asarray(alpha, complex)
        x = np.empty(self.Nv, complex)
        for i, (m, n) in enumerate(self.keys):
            x[i] = np.prod(np.conj(alpha)**np.array(m))*np.prod(alpha**np.array(n))
        return x

    def steady_state(self, x0, T=300.0, rtol=1e-10, atol=1e-12, method='DOP853', blowup=1e4, polish=True):
        out = dict(N=self.N, Nv=self.Nv)
        y0 = np.concatenate([x0.real, x0.imag])
        def rhs(t, y): return self.rhs_real(y)
        def ev(t, y): return blowup - np.max(np.abs(y))
        ev.terminal = True
        t0 = time.time()
        sol = solve_ivp(rhs, (0, T), y0, method=method, rtol=rtol, atol=atol, events=ev)
        out['t_ivp'] = time.time() - t0
        y = sol.y[:, -1]
        out['ivp_status'] = sol.status; out['t_end'] = float(sol.t[-1])
        out['diverged'] = bool(sol.status == 1 or not np.all(np.isfinite(y)))
        out['ivp_residual'] = float(np.max(np.abs(self.rhs_real(y)))) if not out['diverged'] else np.inf
        out['ivp_nfev'] = int(sol.nfev)
        if out['diverged']:
            out['converged'] = False; out['x'] = None; out['t_polish'] = 0.0
            return out
        t0 = time.time()
        if polish:
            ysol, info, ier, msg = fsolve(self.rhs_real, y, full_output=True, xtol=1e-13)
            res = float(np.max(np.abs(self.rhs_real(ysol))))
            if ier == 1 and res < out['ivp_residual']:
                y = ysol
            out['polish_ier'] = int(ier)
        out['t_polish'] = time.time() - t0
        out['residual'] = float(np.max(np.abs(self.rhs_real(y))))
        out['converged'] = out['residual'] < 1e-8
        x = y[:self.Nv] + 1j*y[self.Nv:]
        out['x'] = x
        # hermiticity: <a^dag^m a^n> = conj(<a^dag^n a^m>)
        herm = 0.0
        for (m, n), i in self.idx.items():
            if i == 0: continue
            herm = max(herm, abs(x[i-1] - np.conj(x[self.idx[(n, m)]-1])))
        out['hermiticity'] = float(herm)
        return out

    def jacobian(self, x, h=1e-7):
        y = np.concatenate([x.real, x.imag]); f0 = self.rhs_real(y)
        J = np.empty((len(y), len(y)))
        for i in range(len(y)):
            yp = y.copy(); yp[i] += h
            J[:, i] = (self.rhs_real(yp) - f0)/h
        return J

    def site_stats(self, x, j=0):
        K = self.K
        e = np.array(unit(K, j))
        n = self.moment(x, (tuple(e), tuple(e)))
        n2 = self.moment(x, (tuple(2*e), tuple(2*e)))
        n3 = self.moment(x, (tuple(3*e), tuple(3*e)))
        return dict(n=complex(n), n2=complex(n2), n3=complex(n3),
                    g2=float((n2/n**2).real), g3=float((n3/n**3).real),
                    imag_n=float(abs(n.imag)), imag_n2=float(abs(n2.imag)), imag_n3=float(abs(n3.imag)))


def ring_closure(K, N, kappa, Delta, J, eta, U, ring=True, verbose=False):
    H, jumps = kerr_lattice(K, kappa, Delta, J, eta, U, ring_bonds(K, ring))
    obs = []
    for j in range(K):
        for p in (1, 2, 3):
            e = tuple(int(v) for v in p*np.array(unit(K, j)))
            if 2*p > N: obs.append((e, e))
    return CumulantClosure(K, N, H, jumps, obs_keys=obs, verbose=verbose)


def linear_coherent_alpha(K, kappa, Delta, J, eta, ring=True):
    Om = np.diag(np.asarray(Delta, float)*np.ones(K)).astype(complex)
    for (j, l) in ring_bonds(K, ring):
        Om[j, l] = Om[l, j] = -J
    Gam = kappa/2*np.eye(K) + 1j*Om
    return np.linalg.solve(Gam, np.asarray(eta, complex)*np.ones(K))


if __name__ == "__main__":
    # smoke test: K=3 ring, N=2, U=0.1
    K = 3; P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0)
    for N in (2, 3, 4):
        cc = ring_closure(K, N, U=0.1, verbose=True, **P)
        al = linear_coherent_alpha(K, **P)
        res = cc.steady_state(cc.coherent_state(al))
        st = cc.site_stats(res['x'])
        print(N, res['converged'], res['residual'], res['t_ivp'], res['t_polish'], st['n'].real, st['g2'], st['g3'])
