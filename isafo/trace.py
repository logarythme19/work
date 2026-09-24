"""Simulateur instrumente: meme physique que le noyau, mais avec trajectoires.

Sert au diagnostic et aux figures du par.9.6. Il ne remplace JAMAIS
isafo.evaluate pour decider de la faisabilite: c'est le noyau compile qui fait
foi, et ce module est verifie contre lui.
"""

from __future__ import annotations

import numpy as np

from .params import PLANT, MISSION, NUM
from .evaluate import Controller, _oustaloup_arrays


def trace(ctrl: Controller, h: float = None, plant=None, mission=None,
          stride: int = 1) -> dict:
    """Rejoue la mission en enregistrant t, vo, iL, duty, reference."""
    P = plant if plant is not None else PLANT
    M = mission if mission is not None else MISSION
    h = NUM.h_search if h is None else float(h)
    Tsw = 1.0 / NUM.f_pwm

    zi, pi_, phii, betai, gaini = _oustaloup_arrays(1.0 - ctrl.lam, NUM.Ts)
    zd, pd, phid, betad, gaind = _oustaloup_arrays(ctrl.mu, NUM.Ts)

    n_ctrl = int(round(M.t_end / NUM.Ts))
    iL, vC = M.iL0, M.vC0
    wi = np.zeros(len(pi_)); wd = np.zeros(len(pd))
    xI = 0.0; xF = 0.0
    phiF = np.exp(-ctrl.wf * NUM.Ts)

    dff0 = (P.R + P.rL + P.rd) * M.v0 / (P.R * P.Vin - (P.Ron - P.rd) * M.v0)
    dbuf = np.full(max(NUM.delay_samples, 1), dff0)

    T, VO, IL, D, REF, DT = [], [], [], [], [], []

    def cap(iL, vC, ip):
        ic = (iL - vC / P.R - ip) / (1.0 + P.rC / P.R)
        return ic, vC + P.rC * ic, ic / P.C

    for k in range(n_ctrl):
        t = k * NUM.Ts
        r = M.v0 if t < M.t_step1 else (M.v1 if t < M.t_step2 else M.v2)
        ip = M.ip_amp if (M.ip_t0 <= t < M.ip_t1) else 0.0
        _, vo, _ = cap(iL, vC, ip)

        en = (r - vo) / NUM.V_norm
        eb = (ctrl.b * r - vo) / NUM.V_norm
        ec = (ctrl.c * r - vo) / NUM.V_norm

        y = en
        for i in range(len(pi_)):
            yi = y + (zi[i] - pi_[i]) * wi[i]
            wi[i] = phii[i] * wi[i] + betai[i] * y
            y = yi
        q = gaini * y

        xF = phiF * xF + (1.0 - phiF) * ec
        y = xF
        for i in range(len(pd)):
            yi = y + (zd[i] - pd[i]) * wd[i]
            wd[i] = phid[i] * wd[i] + betad[i] * y
            y = yi
        dterm = gaind * y

        u = ctrl.kp * eb + xI + ctrl.kd * dterm
        dff = (P.R + P.rL + P.rd) * r / (P.R * P.Vin - (P.Ron - P.rd) * r)
        d_tilde = dff + NUM.D_norm * u
        d_sat = min(max(d_tilde, NUM.duty_min), NUM.duty_max)
        xI += NUM.Ts * (ctrl.ki * q + ctrl.Kb * (d_sat - d_tilde) / NUM.D_norm)

        d_app = dbuf[0]
        for s in range(NUM.delay_samples - 1):
            dbuf[s] = dbuf[s + 1]
        dbuf[NUM.delay_samples - 1] = d_sat

        t_local, t_stop = t, t + NUM.Ts
        d_on = d_app * Tsw
        eps = 1e-13
        while t_local < t_stop - eps:
            n_per = np.floor((t_local + eps) / Tsw)
            phase = max(t_local - n_per * Tsw, 0.0)
            if phase < d_on - eps:
                sw_on, t_edge = True, n_per * Tsw + d_on
            else:
                sw_on, t_edge = False, (n_per + 1.0) * Tsw
            if t_edge <= t_local + eps:
                t_edge = t_local + Tsw
            dt = min(t_stop - t_local, t_edge - t_local, h)
            if dt <= eps:
                dt = t_stop - t_local
                if dt <= 0:
                    break

            _, vo_r, _ = cap(iL, vC, ip)
            if sw_on:
                mode = 0
            elif iL > 0.0 or (-P.Vf - vo_r) > 0.0:
                mode = 1
            else:
                mode = 2; iL = 0.0

            def f(i_, v_):
                ic, vo_, dv = cap(i_, v_, ip)
                if mode == 0:
                    di = (P.Vin - vo_ - (P.Ron + P.rL) * i_) / P.L
                elif mode == 1:
                    di = (-P.Vf - (P.rd + P.rL) * i_ - vo_) / P.L
                else:
                    di = 0.0
                return di, dv

            k1 = f(iL, vC); k2 = f(iL + .5*dt*k1[0], vC + .5*dt*k1[1])
            k3 = f(iL + .5*dt*k2[0], vC + .5*dt*k2[1]); k4 = f(iL + dt*k3[0], vC + dt*k3[1])
            iL += dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
            vC += dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
            if mode != 0 and iL < 0.0:
                iL = 0.0

            _, vo_s, _ = cap(iL, vC, ip)
            T.append(t_local + dt); VO.append(vo_s); IL.append(iL)
            D.append(d_app); REF.append(r)
            t_local += dt

    sl = slice(None, None, stride)
    return {"t": np.array(T)[sl], "vo": np.array(VO)[sl], "iL": np.array(IL)[sl],
            "duty": np.array(D)[sl], "ref": np.array(REF)[sl]}
