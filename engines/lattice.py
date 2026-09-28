"""Driven-dissipative Harper-Hofstadter lattice of Kerr resonators.

Gamma = kappa/2 + i Omega with Omega Hermitian (complex Peierls hoppings) is still
normal, so the damped normal modes are orthonormal and the tensor engine applies
unchanged: contractions conserve the mode label, the propagator is diagonal, and the
Kerr term is a dense four-mode vertex.
"""
import numpy as np
import tensor_engine as TE

def harper_hofstadter(Lx, Ly, J, phi, detuning=0.0):
    """Landau gauge: hopping in y picks up exp(i 2 pi phi x)."""
    K = Lx*Ly
    idx = lambda x, y: x*Ly + y
    Om = np.zeros((K, K), complex)
    for x in range(Lx):
        for y in range(Ly):
            i = idx(x, y)
            Om[i, i] = detuning
            if x+1 < Lx:
                j = idx(x+1, y); Om[i, j] += -J; Om[j, i] += -J
            if y+1 < Ly:
                j = idx(x, y+1); ph = np.exp(2j*np.pi*phi*x)
                Om[i, j] += -J*ph; Om[j, i] += -J*np.conj(ph)
    return Om, idx

class Lattice(TE.Chain):
    """Same algebra as Chain, but with an arbitrary Hermitian Omega."""
    def __init__(self, Om, kappa, eta, U):
        K = Om.shape[0]
        self.K, self.U = K, U
        Gam = kappa/2*np.eye(K) + 1j*Om
        eps, Vm = np.linalg.eigh(Om)                 # Hermitian -> orthonormal modes
        self.V = Vm.astype(complex); self.zeta = kappa/2 + 1j*eps
        eta = np.asarray(eta, complex)*np.ones(K)
        self.alpha_site = np.linalg.solve(Gam, eta)
        self.alpha = self.V.conj().T @ self.alpha_site
        self.Hint = self._kerr()
