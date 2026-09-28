"""Verify the general (Wick-ordered) shift relations and the emitter rules."""
import sympy as sp
from math import comb, factorial
N, M, Mc = sp.symbols('N M Mc')   # <a^dag a>=N, <a a>=M, <a^dag a^dag>=Mc ; <a a^dag>=N+1

def nprod(A, B):
    out = {}
    for (m, n), x in A.items():
        for (p, q), y in B.items():
            for k in range(min(n, p)+1):
                key = (m+p-k, n+q-k)
                out[key] = sp.expand(out.get(key, 0) + factorial(k)*comb(n, k)*comb(p, k)*x*y)
    return {k: v for k, v in out.items() if v != 0}
def add(*ops):
    out = {}
    for o in ops:
        for k, v in o.items(): out[k] = sp.expand(out.get(k, 0) + v)
    return {k: v for k, v in out.items() if v != 0}
def sc(c, o): return {k: sp.expand(c*v) for k, v in o.items()}

def dfact2(n):  # number of perfect pairings of 2j objects = (2j)!/(2^j j!)
    return factorial(n)//(2**(n//2)*factorial(n//2)) if n % 2 == 0 else 0
def wick(m, n):
    """:(a^dag)^m a^n: w.r.t. Gaussian state with contractions N, M, Mc, as normal-ordered polynomial"""
    out = {}
    for k in range(min(m, n)+1):                     # mixed pairs a^dag-a, contraction N
        for j in range(0, (m-k)//2+1):               # a^dag-a^dag pairs, contraction Mc
            for l in range(0, (n-k)//2+1):           # a-a pairs, contraction M
                w = comb(m, k)*comb(n, k)*factorial(k) * comb(m-k, 2*j)*dfact2(2*j) * comb(n-k, 2*l)*dfact2(2*l)
                c = (-1)**(k+j+l) * w * N**k * Mc**j * M**l
                key = (m-k-2*j, n-k-2*l)
                out[key] = sp.expand(out.get(key, 0) + c)
    return {kk: v for kk, v in out.items() if v != 0}
a = {(0, 1): 1}; ad = {(1, 0): 1}
def W(m, n): return wick(m, n) if m >= 0 and n >= 0 else {}
ok = True
for m in range(0, 4):
    for n in range(0, 4):
        lhs = nprod(ad, W(m, n)); rhs = add(W(m+1, n), sc(m*Mc, W(m-1, n)), sc(n*N, W(m, n-1)))
        ok &= add(lhs, sc(-1, rhs)) == {}
        lhs = nprod(W(m, n), a); rhs = add(W(m, n+1), sc(m*N, W(m-1, n)), sc(n*M, W(m, n-1)))
        ok &= add(lhs, sc(-1, rhs)) == {}
        lhs = nprod(a, W(m, n)); rhs = add(W(m, n+1), sc(m*(N+1), W(m-1, n)), sc(n*M, W(m, n-1)))
        ok &= add(lhs, sc(-1, rhs)) == {}
        lhs = nprod(W(m, n), ad); rhs = add(W(m+1, n), sc(m*Mc, W(m-1, n)), sc(n*(N+1), W(m, n-1)))
        ok &= add(lhs, sc(-1, rhs)) == {}
print("general boson moves:", ok)

# thermal eigenoperators of L0^dag are :(a^dag)^m a^n: with M=0 (check against paper Eq. 30)
import numpy as np
# emitter: basis I, s-, s+, Xz = s+s- - pe
pe = sp.symbols('p_e')
sm = sp.Matrix([[0, 0], [1, 0]]); spl = sm.T; I = sp.eye(2)
Xz = spl*sm - pe*I
O = {(0, 0): I, (0, 1): sm, (1, 0): spl, (1, 1): Xz}
def expand_basis(Mx):
    # solve Mx = c00 I + c01 s- + c10 s+ + c11 Xz
    c = sp.symbols('c0:4'); expr = c[0]*I + c[1]*sm + c[2]*spl + c[3]*Xz
    sol = sp.solve(list(expr - Mx), c, dict=True)[0]
    return {k: sp.simplify(sol[ci]) for k, ci in zip([(0,0),(0,1),(1,0),(1,1)], c) if sp.simplify(sol[ci]) != 0}
for name, f in [("s+ O", lambda X: spl*X), ("O s-", lambda X: X*sm), ("s- O", lambda X: sm*X), ("O s+", lambda X: X*spl)]:
    print(name, {k: expand_basis(f(v)) for k, v in O.items()})
