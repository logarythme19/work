"""Comparator designs (fixed rules) and decoders of the searched families.

Fixed-rule comparators, all synthesized at the 24 V equilibrium, Ts = 10 us,
same feedforward / saturation / one-sample delay / measurement bus:
  PI    loop shaping: crossover fc = 0.03 fs, phase margin 60 deg including the
        1.5 Ts delay (one computation sample + ZOH), Kb = 1/Ti.
  PID   pole-zero cancellation of the LC pair, phase margin 60 deg including the
        delay and the ESR zero; backward-difference derivative on -vo,
        b = 1, c = 0; Kb = 1/sqrt(Ti Td).
  PIDF  Type-III rule: PID zeros on the LC pair, derivative-filter pole on the
        ESR zero N = 1/(rC C), phase margin 60 deg; b = 1, c = 0; Kb = 1/sqrt(Ti Td).
  LQR   discrete DARE on the ZOH model, Bryson limits of xi0, measured iL, vC.
  LQG   LQR gain + steady-state Kalman predictor on vo only (process noise
        0.01 duty in the Bd direction, measurement noise 0.05 V).
  LQI   discrete DARE on the integral-augmented ZOH model, Bryson limits of
        xi0, measured iL, vC, Kb = wB.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_discrete_are, expm
from scipy.optimize import brentq

from .params import PLANT, SAMP
from .plant import Gvd_coeffs, freq_resp, linearize
from .ctrl import base_vector
from .lqi import lqi, XI0, LB, UB, bryson, DD

TS = SAMP.Ts
DELAY = 1.5 * TS


def _G(w):
    n, d = Gvd_coeffs(24.0)
    return freq_resp(n, d, w)


def pi_rule(fc=0.03 * SAMP.fs, pm=60.0):
    wc = 2 * np.pi * fc
    G = _G([wc])[0]
    phC = np.radians(-180 + pm) - np.angle(G) + wc * DELAY   # required controller phase
    wi = -wc * np.tan(phC)                                     # C = Kp (1 + wi/s): phase = -atan(wi/wc)
    Kp = 1.0 / (abs(G) * abs(1 - 1j * wi / wc))
    Ki = Kp * wi
    p = base_vector(PLANT, SAMP)
    p[0] = 1; p[1] = Kp; p[2] = Ki; p[3] = Ki / Kp
    return p


def pid_rule(pm=60.0, filtered=False):
    n, d = Gvd_coeffs(24.0)
    w0 = np.sqrt(d[2] / d[0]); z0 = (d[1] / d[0]) / (2 * w0)
    Gdc = n[1] / d[2]
    wz = 1.0 / (PLANT.rC * PLANT.C)
    if filtered:
        f = lambda w: 90 - np.degrees(w * DELAY) - pm
        wc = brentq(f, 1.0, np.pi / DELAY)
        Ki = wc / Gdc
    else:
        f = lambda w: 90 + np.degrees(np.arctan(w / wz)) - np.degrees(w * DELAY) - pm
        wc = brentq(f, 1.0, np.pi / DELAY)
        Ki = wc / (Gdc * abs(1 + 1j * wc / wz))
    Kp = Ki * 2 * z0 / w0
    Kd = Ki / w0 ** 2
    Ti, Td = Kp / Ki, Kd / Kp
    Kb = 1.0 / np.sqrt(Ti * Td)
    p = base_vector(PLANT, SAMP)
    if filtered:
        p[0] = 3; p[1:8] = [Kp, Ki, Kd, wz, 1.0, 0.0, Kb]
    else:
        p[0] = 2; p[1:6] = [Kp, Ki, Kd, 1.0, Kb]
    return p, wc


def _zoh(A, B, Ts=TS):
    n = A.shape[0]
    M = np.zeros((n + B.shape[1], n + B.shape[1]))
    M[:n, :n] = A; M[:n, n:] = B
    E = expm(M * Ts)
    return E[:n, :n], E[:n, n:]


def lqr_design():
    lin = linearize(24.0)
    Ad, Bd = _zoh(lin.A, lin.Bd)
    Q, R = bryson(XI0)
    X = solve_discrete_are(Ad, Bd, Q[:2, :2] * TS, R * TS)
    K = np.linalg.solve(R * TS + Bd.T @ X @ Bd, Bd.T @ X @ Ad).ravel()
    return K, Ad, Bd, lin


def lqr_vec():
    K, *_ = lqr_design()
    p = base_vector(PLANT, SAMP)
    p[0] = 4; p[8:10] = K
    return p


def lqg_vec(q_d=0.01, r_v=0.05):
    K, Ad, Bd, lin = lqr_design()
    W = (q_d ** 2) * Bd @ Bd.T
    V = np.array([[r_v ** 2]])
    Cd = lin.Cy
    S = solve_discrete_are(Ad.T, Cd.T, W, V)
    Lf = (S @ Cd.T @ np.linalg.inv(Cd @ S @ Cd.T + V)).ravel()
    p = base_vector(PLANT, SAMP)
    p[0] = 5; p[8:10] = K; p[10:12] = Lf
    p[12:16] = Ad.ravel(); p[16:18] = Bd.ravel(); p[18:20] = Cd.ravel()
    return p


def lqi_vec(xi=XI0):
    lin = linearize(24.0)
    Ad, Bd = _zoh(lin.A, lin.Bd)
    Aa = np.block([[Ad, np.zeros((2, 1))], [-TS * lin.Cy, np.ones((1, 1))]])
    Ba = np.vstack([Bd, [[0.0]]])
    Q, R = bryson(xi)
    X = solve_discrete_are(Aa, Ba, Q * TS, R * TS)
    K = np.linalg.solve(R * TS + Ba.T @ X @ Ba, Ba.T @ X @ Aa).ravel()
    p = base_vector(PLANT, SAMP)
    p[0] = 6; p[1] = 10 ** xi[3]; p[8:11] = K
    return p


# ------------------------------------------------------------ search families
def dec_cmo(x):
    """CMO-LQI-PIDF: Bryson excursions -> CARE -> exact 2-DOF PIDF."""
    th = lqi(x).theta
    p = base_vector(PLANT, SAMP)
    p[0] = 3; p[1:8] = th
    return p


def dec_cmo_sf(x):
    """CMO-LQI, full-state twin: the same continuous gain, measured iL and vC."""
    D = lqi(x)
    p = base_vector(PLANT, SAMP)
    p[0] = 6; p[1] = D.Kb; p[8:11] = D.K
    return p


def _gain_box():
    """Gain ranges spanned by the LQI box (for the direct-search ablations)."""
    rng = np.random.default_rng(0)
    th = np.array([lqi(LB + rng.random(4) * (UB - LB)).theta for _ in range(4000)])
    lo = np.log10(th[:, :3].min(0)); hi = np.log10(th[:, :3].max(0))
    return np.floor(lo * 10) / 10, np.ceil(hi * 10) / 10


GLO, GHI = None, None


def gain_box():
    global GLO, GHI
    if GLO is None:
        GLO, GHI = _gain_box()
    return GLO, GHI


def fam_pidf4():
    lo, hi = gain_box()
    lb = np.r_[lo, LB[3]]; ub = np.r_[hi, UB[3]]
    th0 = lqi(XI0).theta
    x0 = np.r_[np.log10(th0[:3]), XI0[3]]
    def dec(x):
        p = base_vector(PLANT, SAMP)
        p[0] = 3
        p[1:8] = [10 ** x[0], 10 ** x[1], 10 ** x[2], 1 / (PLANT.rC * PLANT.C), 1.0, 0.0, 10 ** x[3]]
        return p
    return dec, lb, ub, x0


def fam_pidf7():
    lo, hi = gain_box()
    lb = np.r_[lo, 3.0, 0.0, 0.0, LB[3]]
    ub = np.r_[hi, np.log10(np.pi / TS), 1.0, 1.0, UB[3]]
    th0 = lqi(XI0).theta
    x0 = np.r_[np.log10(th0[:4]), 1.0, 0.0, XI0[3]]
    def dec(x):
        p = base_vector(PLANT, SAMP)
        p[0] = 3
        p[1:8] = [10 ** x[0], 10 ** x[1], 10 ** x[2], 10 ** x[3], x[4], x[5], 10 ** x[6]]
        return p
    return dec, lb, ub, x0


def fam_pid4():
    lo, hi = gain_box()
    lb = np.r_[lo, LB[3]]; ub = np.r_[hi, UB[3]]
    p0, _ = pid_rule()
    x0 = np.clip(np.log10([p0[1], p0[2], p0[3], p0[5]]), lb, ub)
    def dec(x):
        p = base_vector(PLANT, SAMP)
        p[0] = 2; p[1:6] = [10 ** x[0], 10 ** x[1], 10 ** x[2], 1.0, 10 ** x[3]]
        return p
    return dec, lb, ub, x0


def fam_pi3():
    lo, hi = gain_box()
    lb = np.r_[lo[:2], LB[3]]; ub = np.r_[hi[:2], UB[3]]
    p0 = pi_rule()
    x0 = np.clip(np.log10([p0[1], p0[2], p0[3]]), lb, ub)
    def dec(x):
        p = base_vector(PLANT, SAMP)
        p[0] = 1; p[1:4] = [10 ** x[0], 10 ** x[1], 10 ** x[2]]
        return p
    return dec, lb, ub, x0


FAMILIES = {
    'CMO-LQI-PIDF': lambda: (dec_cmo, LB, UB, XI0),
    'CMO-LQI-SF': lambda: (dec_cmo_sf, LB, UB, XI0),
    'PIDF-direct': fam_pidf4,
    'PIDF7-direct': fam_pidf7,
    'PID-direct': fam_pid4,
    'PI-direct': fam_pi3,
}


def fixed_comparators():
    return {
        'PI': pi_rule(),
        'PID': pid_rule(filtered=False)[0],
        'PIDF': pid_rule(filtered=True)[0],
        'LQR': lqr_vec(),
        'LQG': lqg_vec(),
        'LQI': lqi_vec(),
    }
