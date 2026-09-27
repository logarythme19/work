"""Trajectory-level test of the LQI <-> 2-DOF PIDF equivalence (Proposition 3).

Continuous-time controllers on the nonlinear averaged plant (Eq. (6)) with
duty saturation, back-calculation anti-windup and the iL >= 0 guard:
  law 0 (LQI, measured iL and vC):
      d~ = dff(r) - Kix (iL - r/R) - Kvx (vC - r) - Kz z,
      z' = (r - vo) + Kb/|Kz| (d - d~)
  law 1 (2-DOF PIDF, measured vo only, b = 1, c = 0, N = 1/(rC C)):
      d~ = dff(r) + Kp (r - vo) + xI - Kd N (vo + xF),
      xI' = Ki (r - vo) + Kb (d - d~),   xF' = N(-vo - xF),   xF(0) = -vC(0)
Proposition 3 predicts identical trajectories whenever ip = 0 and the output
network (R, rC, C) equals its nominal value, whatever the inductor branch,
saturation or anti-windup activity.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _f(x, law, r, ip, pl, Rn, K, th, Kb, dmin, dmax):
    Vin, L, C, R, rL, rC, Ron, rd, Vf = pl[0], pl[1], pl[2], pl[3], pl[4], pl[5], pl[6], pl[7], pl[8]
    iL, vC, xi, xf = x[0], x[1], x[2], x[3]
    ic = (iL - vC / R - ip) / (1.0 + rC / R)
    vo = vC + rC * ic
    ILn = r / Rn
    dff = (r + (rL + rd) * ILn + Vf) / (Vin - (Ron - rd) * ILn + Vf)
    if law == 0:
        dt = dff - K[0] * (iL - r / Rn) - K[1] * (vC - r) - K[2] * xi
    else:
        dt = dff + th[0] * (r - vo) + xi - th[2] * th[3] * (vo + xf)
    d = min(max(dt, dmin), dmax)
    diL = (d * Vin - (rL + d * Ron + (1 - d) * rd) * iL - (1 - d) * Vf - vo) / L
    if iL <= 0.0 and diL < 0.0:
        diL = 0.0
    out = np.empty(4)
    out[0] = diL
    out[1] = ic / C
    if law == 0:
        out[2] = (r - vo) + Kb / abs(K[2]) * (d - dt)
        out[3] = 0.0
    else:
        out[2] = th[1] * (r - vo) + Kb * (d - dt)
        out[3] = th[3] * (-vo - xf)
    return out, vo, d, dt


@njit(cache=True)
def sim_cont(law, pl, Rn, K, th, Kb, tr, lv, ld_on, ld_off, ld_amp, iL0, vC0, T, h, dmin, dmax, nrec):
    n = int(round(T / h))
    x = np.array([iL0, vC0, 0.0, -vC0])
    every = max(n // nrec, 1)
    m = n // every + 1
    rec = np.zeros((m, 5))
    j = 0
    for k in range(n + 1):
        t = k * h
        r = lv[0]
        for q in range(tr.shape[0]):
            if t >= tr[q]:
                r = lv[q + 1]
        ip = ld_amp if (ld_on <= t < ld_off) else 0.0
        if k % every == 0 and j < m:
            _, vo, d, dt = _f(x, law, r, ip, pl, Rn, K, th, Kb, dmin, dmax)
            rec[j, 0] = t; rec[j, 1] = vo; rec[j, 2] = x[0]; rec[j, 3] = d; rec[j, 4] = x[1]
            j += 1
        if k == n:
            break
        k1, _, _, _ = _f(x, law, r, ip, pl, Rn, K, th, Kb, dmin, dmax)
        k2, _, _, _ = _f(x + 0.5 * h * k1, law, r, ip, pl, Rn, K, th, Kb, dmin, dmax)
        k3, _, _, _ = _f(x + 0.5 * h * k2, law, r, ip, pl, Rn, K, th, Kb, dmin, dmax)
        k4, _, _, _ = _f(x + h * k3, law, r, ip, pl, Rn, K, th, Kb, dmin, dmax)
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if x[0] < 0.0:
            x[0] = 0.0
    return rec[:j]
