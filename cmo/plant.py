"""Formula-first non-ideal Buck model: Eqs. (1)-(15) of the article.

Everything here is closed form; the numerical instantiations reported in the
paper (Tables 2 and 3, Eq. (16)) are produced by scripts/cmo_formula_first.py.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .params import PLANT, SAMP


def gamma(P=PLANT):
    """gamma = R / (R + rC), Eq. (8)."""
    return P.R / (P.R + P.rC)


def equilibrium(Vo, Ip=0.0, P=PLANT):
    """Loss-aware equilibrium duty, Eq. (7). Returns (IL, d0)."""
    IL = Vo / P.R + Ip
    d0 = (Vo + (P.rL + P.rd) * IL + P.Vf) / (P.Vin - (P.Ron - P.rd) * IL + P.Vf)
    return IL, d0


def dff(r, P=PLANT):
    """Nominal feedforward used online by every controller (ip = 0)."""
    return equilibrium(r, 0.0, P)[1]


@dataclass
class Linear:
    A: np.ndarray
    Bd: np.ndarray
    Bg: np.ndarray
    Bp: np.ndarray
    Cy: np.ndarray
    Dp: float
    d0: float
    IL0: float
    rho0: float
    Ed: float
    gamma: float


def linearize(Vo, Ip=0.0, P=PLANT) -> Linear:
    """Small-signal model, Eqs. (8)-(10)."""
    IL, d0 = equilibrium(Vo, Ip, P)
    g = gamma(P)
    rho0 = P.rL + P.rd + d0 * (P.Ron - P.rd)
    Ed = P.Vin + P.Vf - (P.Ron - P.rd) * IL
    A = np.array([[-(rho0 + g * P.rC) / P.L, -g / P.L],
                  [g / P.C, -g / (P.R * P.C)]])
    Bd = np.array([[Ed / P.L], [0.0]])
    Bg = np.array([[d0 / P.L], [0.0]])
    Bp = np.array([[g * P.rC / P.L], [-g / P.C]])
    Cy = np.array([[g * P.rC, g]])
    return Linear(A, Bd, Bg, Bp, Cy, -g * P.rC, d0, IL, rho0, Ed, g)


def Gvd_coeffs(Vo, Ip=0.0, P=PLANT):
    """Closed form of Eqs. (11)-(12): numerator and denominator (descending s)."""
    lin = linearize(Vo, Ip, P)
    rho0 = lin.rho0
    num = lin.Ed * P.R * np.array([P.rC * P.C, 1.0])
    den = np.array([P.L * P.C * (P.R + P.rC),
                    P.L + P.C * ((P.R + P.rC) * rho0 + P.R * P.rC),
                    P.R + rho0])
    return num, den


def Gvg_coeffs(Vo, Ip=0.0, P=PLANT):
    """Eq. (13)."""
    lin = linearize(Vo, Ip, P)
    num = lin.d0 * P.R * np.array([P.rC * P.C, 1.0])
    return num, Gvd_coeffs(Vo, Ip, P)[1]


def ccm_screen(Vo, P=PLANT, fs=SAMP.fs):
    """Stationary CCM/DCM screen, Eq. (15). Returns a dict."""
    IL, d = equilibrium(Vo, 0.0, P)
    dI = Vo * (1 - d) / (P.L * fs)
    Icrit = dI / 2
    Rcrit = 2 * P.L * fs / (1 - d)
    return dict(Vo=Vo, dff=d, dIpp=dI, Icrit=Icrit, margin=IL - Icrit,
                margin_pct=100 * (IL - Icrit) / IL, Rcrit=Rcrit,
                Rcrit_over_R=Rcrit / P.R)


def output_ripple_pp(Vo, P=PLANT, fs=SAMP.fs):
    """Peak-to-peak output-voltage ripple at the carrier frequency (CCM).

    ESR term plus capacitive term: dVo = dI (rC + 1/(8 C fs)). Used by the
    switching-ripple constraint g37 (see cmo/ripple.py).
    """
    IL, d = equilibrium(Vo, 0.0, P)
    dI = Vo * (1 - d) / (P.L * fs)
    return dI * np.hypot(P.rC, 1.0 / (8 * P.C * fs))


def freq_resp(num, den, w):
    s = 1j * np.asarray(w)
    return np.polyval(num, s) / np.polyval(den, s)
