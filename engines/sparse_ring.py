"""Sparse-matrix form of the momentum-space ring engine (multimode.Ring).

The recursion X_N = V G_0 X_{N-1} of multimode.Ring is linear in X, so the map V G_0 is a matrix on
the space of labels.  multimode.Ring rebuilds that map symbolically at every order, for every
observable and for every value of U.  Here it is built once:

  * the reachable labels are enumerated layer by layer from the observable, keeping at step N only
    labels of degree <= drop*(Nmax-N) (a Kerr vertex lowers the degree by at most two, so higher
    labels cannot return to the steady state) and of zero net momentum (the vertex conserves
    momentum and only zero-momentum labels can reach the (0,0) label);
  * the column of V G_0 for each label is computed once with the dictionary algebra of multimode
    (at U = 1) and stored as a scipy.sparse matrix;
  * the series is then Nmax sparse matrix-vector products, and since the vertex is proportional to U
    the N-th order term is c_N U^N: one run gives the partial sums at every U.

    R = SparseRing(K=3, kappa=1, Delta=-1, J=0.4, eta=1)       # no U
    b, bd = R.site_ops(0)
    c = R.coefficients(MM.mul(bd, b), Nmax=10)                # c[N] = <b^dag b>^(N) at U = 1
    S = R.partial_sums(c, U)                                  # partial sums S_N(U) = sum_{n<=N} c_n U^n

Validated against multimode.Ring / ring_trunc.series_pruned to machine precision (see check_sparse_ring.py).
"""
import numpy as np
import scipy.sparse as sps
import multimode as MM

def degree(k):
    return sum(c for _, c in k[0]) + sum(c for _, c in k[1])

class SparseRing:
    def __init__(self, K, kappa, Delta, J, eta, eps=None, drop=2, momentum_filter=True):
        self.ring = MM.Ring(K, kappa, Delta, J, eta, 1.0, eps=eps)      # U = 1: the vertex is linear in U
        self.K, self.drop, self.momentum_filter = K, drop, momentum_filter
        self._col = {}                                                   # label -> {label': coefficient}
        self.alpha = self.ring.alpha

    def site_ops(self, j=0):
        return self.ring.site_ops(j)

    def net_momentum(self, k):
        return (sum(q*c for q, c in k[1]) - sum(q*c for q, c in k[0])) % self.K

    def column(self, k):
        """the image of the single label k under V G_0 (computed once, cached)"""
        if k not in self._col:
            self._col[k] = self.ring.V(self.ring.G0({k: 1.0}))
        return self._col[k]

    def labels(self, O, Nmax):
        """layered label sets: layer[N] = labels that can be present after vertex N"""
        keep = lambda k, N: degree(k) <= self.drop*(Nmax-N) and (not self.momentum_filter or self.net_momentum(k) == 0)
        layer = [{k for k in O if keep(k, 0)}]
        for N in range(1, Nmax+1):
            new = set()
            for k in layer[-1]:
                new.update(kk for kk in self.column(k) if keep(kk, N))
            layer.append(new)
        return layer

    def matrix(self, O, Nmax):
        """index of all reachable labels and the sparse matrix M of V G_0 restricted to them"""
        layer = self.labels(O, Nmax)
        idx = {}
        for L in layer:
            for k in L:
                if k not in idx: idx[k] = len(idx)
        rows, cols, vals = [], [], []
        for k, i in idx.items():
            for kk, v in self.column(k).items():
                if kk in idx:
                    rows.append(idx[kk]); cols.append(i); vals.append(v)
        M = sps.csc_matrix((np.array(vals, complex), (rows, cols)), shape=(len(idx), len(idx)))
        masks = []
        for L in layer:
            m = np.zeros(len(idx), bool)
            for k in L: m[idx[k]] = True
            masks.append(m)
        return idx, M, masks

    def coefficients(self, O, Nmax, return_counts=False):
        """c[N] = coefficient of U^N in <O>, N = 0..Nmax, with the same pruning as ring_trunc.series_pruned"""
        idx, M, masks = self.matrix(O, Nmax)
        i00 = idx.get(((), ()), None)
        x = np.zeros(len(idx), complex)
        for k, v in O.items():
            if k in idx: x[idx[k]] = v
        x *= masks[0]
        c = [x[i00] if i00 is not None else 0.0]
        for N in range(1, Nmax+1):
            x = M @ x
            x *= masks[N]
            c.append(x[i00] if i00 is not None else 0.0)
        c = np.array(c)
        return (c, [int(m.sum()) for m in masks]) if return_counts else c

    @staticmethod
    def partial_sums(c, U):
        """S_N(U) = sum_{n<=N} c_n U^n for a scalar U or an array of U (shape (len(U), Nmax+1))"""
        U = np.asarray(U, float)
        pw = U[..., None]**np.arange(len(c))
        return np.cumsum(pw*c, axis=-1)
