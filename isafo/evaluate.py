"""Evaluateur de mission: theta FO-PIDF -> (J, g[6], marge, diagnostics).

Le cout et les six contraintes sont ceux du par.1 du cadrage, sans
reponderation ni assouplissement:

    J = int|r-v|/max(r,1) dt / 0.324
      + 0.01 * int iL^2 dt / (0.324 * 20^2)
      + 0.001 * mean|delta d| / 0.96

    g1 = Imax/20 - 1          g4 = |e1|/0.03  - 1
    g2 = (Vmax-46)/1.5 - 1    g5 = |e46|/0.10 - 1
    g3 = |e24|/0.05 - 1       g6 = (0.8 - Vmin)/0.8

La marge est m = -max_i g_i. C'est un min de max, donc NON DIFFERENTIABLE
(par.3.4): tout raffinement local doit en tenir compte.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .params import PLANT, MISSION, NUM, COST
from .oustaloup import oustaloup_zpk
from ._kernel import simulate


@dataclass
class Controller:
    """Correcteur FO-PIDF en coordonnees normalisees.

    theta = [log10(kp), log10(ki), log10(kd), lambda, mu, log10(wf), b], c = 0.
    Kb (anti-windup) provient du vecteur de synthese, pas de theta.
    """

    kp: float
    ki: float
    kd: float
    lam: float
    mu: float
    wf: float
    b: float
    Kb: float
    c: float = 0.0

    @classmethod
    def from_theta(cls, theta, Kb: float) -> "Controller":
        th = np.asarray(theta, dtype=float)
        if th.shape[0] != 7:
            raise ValueError(f"theta doit avoir 7 composantes, recu {th.shape[0]}")
        return cls(kp=10.0 ** th[0], ki=10.0 ** th[1], kd=10.0 ** th[2],
                   lam=float(th[3]), mu=float(th[4]), wf=10.0 ** th[5],
                   b=float(th[6]), Kb=float(Kb))


@dataclass
class Result:
    J: float
    g: np.ndarray
    margin: float
    feasible: bool
    violation: float
    Imax: float
    Vmax: float
    Vmin: float
    e24: float
    e1: float
    e46: float
    dcm_fraction: float
    diverged: bool
    h: float

    def active_constraint(self) -> int:
        """Indice (1..6) de la contrainte active, c.-a-d. la plus serree."""
        return int(np.argmax(self.g)) + 1


_OUST_CACHE: dict = {}


def _oustaloup_arrays(gamma: float, h_ctrl: float):
    """Tableaux discrets d'une cascade d'Oustaloup, memoises par (gamma, Ts)."""
    key = (round(float(gamma), 12), round(float(h_ctrl), 15))
    hit = _OUST_CACHE.get(key)
    if hit is not None:
        return hit
    z, p, gain = oustaloup_zpk(gamma)
    phi = np.exp(-p * h_ctrl)
    beta = (1.0 - phi) / p
    out = (np.ascontiguousarray(z), np.ascontiguousarray(p),
           np.ascontiguousarray(phi), np.ascontiguousarray(beta), float(gain))
    _OUST_CACHE[key] = out
    return out


def evaluate(ctrl: Controller,
             h: Optional[float] = None,
             plant=None,
             mission=None) -> Result:
    """Execute la mission complete et renvoie cout, contraintes et marge.

    `h` vaut par defaut le pas de recherche (1 us). Le replay fin (0.2 us)
    s'obtient en passant h = NUM.h_confirm; aucun candidat n'est declare
    faisable sans ce replay (par.5).
    """
    P = plant if plant is not None else PLANT
    M = mission if mission is not None else MISSION
    h = NUM.h_search if h is None else float(h)

    # Voie integrale: I^lambda = (1/s) * s^(1-lambda) -> ordre 1-lambda.
    zi, pi_, phii, betai, gaini = _oustaloup_arrays(1.0 - ctrl.lam, NUM.Ts)
    # Voie derivee: s^mu, realisation derivee uniforme (active meme a mu = 1).
    zd, pd, phid, betad, gaind = _oustaloup_arrays(ctrl.mu, NUM.Ts)

    Tsw = 1.0 / NUM.f_pwm

    out = simulate(
        P.Vin, P.L, P.C, P.R, P.rL, P.rC, P.Ron, P.rd, P.Vf,
        M.v0, M.v1, M.v2, M.t_step1, M.t_step2, M.t_end,
        M.ip_amp, M.ip_t0, M.ip_t1, M.iL0, M.vC0,
        NUM.Ts, h, Tsw, NUM.duty_min, NUM.duty_max, NUM.delay_samples,
        NUM.V_norm, NUM.D_norm,
        PLANT.R, PLANT.rL, PLANT.rd, PLANT.Ron, PLANT.Vin,
        ctrl.kp, ctrl.ki, ctrl.kd, ctrl.wf, ctrl.b, ctrl.c, ctrl.Kb,
        zi, pi_, phii, betai, gaini,
        zd, pd, phid, betad, gaind,
        M.win_24[0], M.win_24[1], M.win_1[0], M.win_1[1],
        M.win_46[0], M.win_46[1],
    )

    (int_err, int_iL2, sum_dd, n_duty,
     Imax, Vmax, Vmin, e24, e1, e46,
     n_dcm, n_sub, diverged) = out

    T = M.t_end

    if diverged or n_duty == 0:
        g = np.full(COST.n_constraints, 1e6)
        return Result(J=1e6, g=g, margin=-1e6, feasible=False, violation=6e6,
                      Imax=float(Imax), Vmax=float(Vmax), Vmin=float(Vmin),
                      e24=float(e24), e1=float(e1), e46=float(e46),
                      dcm_fraction=float(n_dcm) / max(n_sub, 1),
                      diverged=True, h=h)

    J = (int_err / T
         + COST.w_current * int_iL2 / (T * COST.i_ref ** 2)
         + COST.w_duty * (sum_dd / n_duty) / COST.dduty_ref)

    g = np.empty(COST.n_constraints)
    g[0] = Imax / COST.imax_lim - 1.0
    g[1] = (Vmax - COST.vmax_ref) / COST.vmax_band - 1.0
    g[2] = abs(e24) / COST.e24_tol - 1.0
    g[3] = abs(e1) / COST.e1_tol - 1.0
    g[4] = abs(e46) / COST.e46_tol - 1.0
    g[5] = (COST.vmin_lim - Vmin) / COST.vmin_lim

    violation = float(np.sum(np.maximum(g, 0.0)))
    margin = float(-np.max(g))

    return Result(J=float(J), g=g, margin=margin,
                  feasible=bool(np.all(g <= 0.0)), violation=violation,
                  Imax=float(Imax), Vmax=float(Vmax), Vmin=float(Vmin),
                  e24=float(e24), e1=float(e1), e46=float(e46),
                  dcm_fraction=float(n_dcm) / max(n_sub, 1),
                  diverged=False, h=h)
