"""Approximation d'Oustaloup de s^gamma et sa realisation discrete exacte.

Bande et ordre imposes par le par.1 du cadrage: demi-ordre 3 -> 2*3+1 = 7
sections, bande 1..1e5 rad/s.

    s^gamma  ~=  K * prod_{k=-N}^{N} (s + w'_k) / (s + w_k)
    K        =  wh^gamma
    w'_k     =  wb * (wh/wb)^((k+N+0.5(1-gamma))/(2N+1))
    w_k      =  wb * (wh/wb)^((k+N+0.5(1+gamma))/(2N+1))

La formule reste bien definie aux deux bouts: gamma = 0 donne exactement 1, et
gamma = 1 donne wh (s+wb)/(s+wh), c'est-a-dire une derivee a bande limitee.
C'est precisement la "realisation derivee uniforme" exigee par le cadrage, qui
reste donc active pour mu = 1 au lieu de basculer sur une derivee exacte.

Chaque section (s+z)/(s+p) s'ecrit 1 + (z-p)/(s+p) et se discretise EXACTEMENT
par bloqueur d'ordre zero, ce qui evite toute integration numerique du
correcteur et rend le filtre inconditionnellement stable quel que soit le pas.
"""

from __future__ import annotations

import numpy as np

from .params import NUM


def oustaloup_zpk(gamma: float,
                  N: int = None,
                  wb: float = None,
                  wh: float = None):
    """Zeros, poles et gain de l'approximation d'Oustaloup de s^gamma.

    gamma doit appartenir a [0, 1]. Les ordres superieurs a 1 doivent etre
    factorises par l'appelant (s^1.3 = s * s^0.3), ce que la chaine FO-PIDF
    n'utilise pas: mu <= 1 et 1-lambda <= 1 par construction.
    """
    if not (0.0 - 1e-12 <= gamma <= 1.0 + 1e-12):
        raise ValueError(f"gamma hors de [0,1]: {gamma}")
    gamma = float(np.clip(gamma, 0.0, 1.0))

    N = NUM.oustaloup_half_order if N is None else N
    wb = NUM.oustaloup_wb if wb is None else wb
    wh = NUM.oustaloup_wh if wh is None else wh

    n_sec = 2 * N + 1
    k = np.arange(-N, N + 1, dtype=float)
    ratio = wh / wb
    zeros = wb * ratio ** ((k + N + 0.5 * (1.0 - gamma)) / n_sec)
    poles = wb * ratio ** ((k + N + 0.5 * (1.0 + gamma)) / n_sec)
    gain = wh ** gamma
    return zeros, poles, gain


class FractionalOperator:
    """Realisation discrete a pas fixe de s^gamma (bloqueur d'ordre zero).

    Etat interne: un scalaire par section. La mise a jour est exacte, donc
    independante de la raideur du filtre (poles jusqu'a 1e5 rad/s pour
    Ts = 10 us).
    """

    __slots__ = ("gamma", "zeros", "poles", "gain", "h", "_phi", "_beta", "w")

    def __init__(self, gamma: float, h: float,
                 N: int = None, wb: float = None, wh: float = None):
        self.gamma = float(gamma)
        self.h = float(h)
        self.zeros, self.poles, self.gain = oustaloup_zpk(gamma, N, wb, wh)

        p = self.poles
        self._phi = np.exp(-p * self.h)                 # e^{-p h}
        self._beta = (1.0 - self._phi) / p              # integrale du BOZ
        self.w = np.zeros_like(p)

    def reset(self) -> None:
        self.w[:] = 0.0

    def step(self, u: float) -> float:
        """Avance d'un pas h et renvoie la sortie de s^gamma appliquee a u."""
        y = u
        w = self.w
        z, p = self.zeros, self.poles
        phi, beta = self._phi, self._beta
        for i in range(len(p)):
            yi = y + (z[i] - p[i]) * w[i]
            w[i] = phi[i] * w[i] + beta[i] * y
            y = yi
        return self.gain * y

    def frequency_response(self, omega: np.ndarray) -> np.ndarray:
        """Reponse frequentielle continue de l'approximation (pour audit)."""
        s = 1j * np.asarray(omega, dtype=complex)
        num = np.ones_like(s)
        den = np.ones_like(s)
        for z, p in zip(self.zeros, self.poles):
            num = num * (s + z)
            den = den * (s + p)
        return self.gain * num / den
