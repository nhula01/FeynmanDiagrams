"""Inverse of the Wick construction: express (c^dag)^m c^n in the eigenoperator basis
O^{(N,M)}_{pq}, and the displaced version (a^dag)^m a^n with a = b + alpha."""
import numpy as np
from math import factorial, comb
from engine import destroy, steady, cav, Model, add, scale, mul

def wick_coeffs(m, n, N, M, sign=-1):
    """sign=-1: O_mn in terms of monomials.  sign=+1: monomials in terms of O_pq."""
    out = {}
    for p in range(m//2+1):
        for q in range(n//2+1):
            for k in range(min(m-2*p, n-2*q)+1):
                c = (factorial(m)*factorial(n)*(sign*np.conj(M))**p*(sign*M)**q*(sign*N)**k
                     /(2**(p+q)*factorial(p)*factorial(q)*factorial(k)
                       *factorial(m-2*p-k)*factorial(n-2*q-k)))
                out[(m-2*p-k, n-2*q-k)] = out.get((m-2*p-k, n-2*q-k), 0) + c
    return out

if __name__ == "__main__":
    Nc = 70
    a = destroy(Nc); ad = a.conj().T
    omega, gam, rb, th, kp, km = 0.7, 1.0, 0.45, 0.3, 0.15, 0.35
    L = np.cosh(rb)*a + np.exp(1j*th)*np.sinh(rb)*ad
    H = omega*ad@a; cops = [np.sqrt(gam)*L, np.sqrt(kp)*ad, np.sqrt(km)*a]
    rho = steady(H, cops)
    N = np.trace(ad@a@rho); M = np.trace(a@a@rho)
    mp = lambda A, k: np.linalg.matrix_power(A, k)
    def O(m, n):
        out = np.zeros((Nc, Nc), complex)
        for (p, q), c in wick_coeffs(m, n, N, M, -1).items():
            out += c*mp(ad, p)@mp(a, q)
        return out
    worst = 0; K = 12
    for m in range(4):
        for n in range(4):
            lhs = mp(ad, m)@mp(a, n)
            rhs = sum(c*O(p, q) for (p, q), c in wick_coeffs(m, n, N, M, +1).items())
            worst = max(worst, np.max(np.abs((lhs-rhs)[:K, :K])))
            c00 = wick_coeffs(m, n, N, M, +1).get((0, 0), 0)
            assert abs(c00 - np.trace(lhs@rho)) < 1e-8, (m, n, c00, np.trace(lhs@rho))
    print("inversion residual:", worst, " and O_00 coefficient = steady moment: ok")
    print("examples: <a^dag a> ->", wick_coeffs(1,1,N,M,+1), "\n  (a^dag)^2a^2 ->",
          {k: np.round(v,4) for k, v in wick_coeffs(2,2,N,M,+1).items()})

def displaced_moment(mod, m, n, alpha, Nmax):
    """<(a^dag)^m a^n> to order Nmax, from series for O_{jk}=(b^dag)^j b^k (vacuum: N=M=0)."""
    tot = np.zeros(Nmax+1, complex)
    for j in range(m+1):
        for k in range(n+1):
            c = comb(m, j)*comb(n, k)*np.conj(alpha)**(m-j)*alpha**(n-k)
            tot += c*mod.series(cav(j, k, np.eye(1)), Nmax)
    return tot
