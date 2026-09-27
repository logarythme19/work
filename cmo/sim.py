"""Averaged (search) and switched (validation) simulators, numba-compiled.

Averaged model  : Eq. (6) with the hybrid guard iL >= 0 (averaged DCM), RK4.
Switched model  : replica of the reference .slx power stage (Fig. 2):
  - high-side switch: resistance Ron when on, resistive snubber Rs always;
  - freewheel diode: affine law Vf + rd*iD, no snubber, blocks reverse current
    (discontinuous conduction arises naturally; the stiff snubber-leakage
    branch is solved quasi-statically, time constant L/Rs = 1 ns);
  - PWM: sawtooth carrier at fs, q = 1 while phase < d*Tsw; edges are located
    exactly inside the integration step;
  - controller sampled at Ts with one-sample delay (S-function semantics);
  - energy boundary of ENERGY_ACCOUNTING in the .slx:
      Pin = Vin*i_in, Pu = vo^2/R + vo*ip,
      Ploss = rL iL^2 + rC iC^2 + vSW iSW + vD iD   (snubber inside vSW iSW)
    integrated by the same RK4 stages as the states, so the residual
    Ein - Eu - Eloss - dEstored measures integration consistency only.
Both simulators return the Ts-sampled traces and scalar diagnostics.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from .ctrl import ctrl_init, ctrl_step, NS


@njit(cache=True)
def _ref(t, tr, lv):
    r = lv[0]
    for k in range(tr.shape[0]):
        if t >= tr[k]:
            r = lv[k + 1]
    return r


@njit(cache=True)
def _cap(iL, vC, ip, R, rC):
    ic = (iL - vC / R - ip) / (1.0 + rC / R)
    return ic, vC + rC * ic


# ------------------------------------------------------------------ averaged
@njit(cache=True)
def _avg_f(iL, vC, d, ip, pl):
    Vin, L, C, R, rL, rC, Ron, rd, Vf = pl[0], pl[1], pl[2], pl[3], pl[4], pl[5], pl[6], pl[7], pl[8]
    ic, vo = _cap(iL, vC, ip, R, rC)
    diL = (d * Vin - (rL + d * Ron + (1 - d) * rd) * iL - (1 - d) * Vf - vo) / L
    if iL <= 0.0 and diL < 0.0:
        diL = 0.0
    pin = d * Vin * iL
    pu = vo * vo / R + vo * ip
    pl_ = (rL + d * Ron + (1 - d) * rd) * iL * iL + (1 - d) * Vf * iL + rC * ic * ic
    return diL, ic / C, pin, pu, pl_


@njit(cache=True)
def sim_avg(pl, p, tr, lv, ld_on, ld_off, ld_amp, iL0, vC0, T, nsub):
    """Averaged simulation. Returns traces at Ts and energies."""
    Ts = p[32]
    n = int(round(T / Ts))
    h = Ts / nsub
    s = np.zeros(NS)
    iL = iL0; vC = vC0
    r0 = _ref(0.0, tr, lv)
    ic, vo = _cap(iL, vC, 0.0, pl[3], pl[5])
    ctrl_init(p, s, r0, vo, iL, vC)
    tv = np.zeros(n + 1); rv = np.zeros(n + 1); vv = np.zeros(n + 1)
    iv = np.zeros(n + 1); dv = np.zeros(n + 1); dpre = np.zeros(n + 1)
    Ein = 0.0; Eu = 0.0; El = 0.0
    nclip = 0; nsat = 0; ipk = iL
    for k in range(n):
        t = k * Ts
        r = _ref(t, tr, lv)
        ip = ld_amp if (ld_on <= t < ld_off) else 0.0
        ic, vo = _cap(iL, vC, ip, pl[3], pl[5])
        dnow, dsat, dt_ = ctrl_step(p, s, r, vo, iL, ic)
        tv[k] = t; rv[k] = r; vv[k] = vo; iv[k] = iL; dv[k] = dnow; dpre[k] = dt_
        if dsat != dt_:
            nsat += 1
        for j in range(nsub):
            a1, b1, p1, u1, l1 = _avg_f(iL, vC, dnow, ip, pl)
            a2, b2, p2, u2, l2 = _avg_f(iL + 0.5 * h * a1, vC + 0.5 * h * b1, dnow, ip, pl)
            a3, b3, p3, u3, l3 = _avg_f(iL + 0.5 * h * a2, vC + 0.5 * h * b2, dnow, ip, pl)
            a4, b4, p4, u4, l4 = _avg_f(iL + h * a3, vC + h * b3, dnow, ip, pl)
            iL += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            vC += h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
            Ein += h / 6 * (p1 + 2 * p2 + 2 * p3 + p4)
            Eu += h / 6 * (u1 + 2 * u2 + 2 * u3 + u4)
            El += h / 6 * (l1 + 2 * l2 + 2 * l3 + l4)
            if iL < 0.0:
                iL = 0.0
                nclip += 1
            if iL > ipk:
                ipk = iL
        if not (abs(vC) < 1e4 and abs(iL) < 1e4):
            return tv, rv, vv, iv, dv, dpre, np.array([np.inf, np.inf, np.inf, 1.0, 1.0, 1e9, 1.0])
    t = n * Ts
    r = _ref(t, tr, lv)
    ip = ld_amp if (ld_on <= t < ld_off) else 0.0
    ic, vo = _cap(iL, vC, ip, pl[3], pl[5])
    tv[n] = t; rv[n] = r; vv[n] = vo; iv[n] = iL; dv[n] = s[4]; dpre[n] = s[5]
    diag = np.array([Ein, Eu, El, nclip / (n * nsub), nsat / n, ipk, 0.0])
    return tv, rv, vv, iv, dv, dpre, diag


# ------------------------------------------------------------------ switched
@njit(cache=True)
def _sw_f(iL, vC, q, ip, pl, fm):
    """Returns diL, dvC, pin, pu, ploss, mode (0 on, 1 diode, 2 DCM)."""
    Vin, L, C, R, rL, rC, Ron, rd, Vf, Rs = pl[0], pl[1], pl[2], pl[3], pl[4], pl[5], pl[6], pl[7], pl[8], pl[9]
    ic, vo = _cap(iL, vC, ip, R, rC)
    Rb = pl[10]
    ith = (Vin + Vf) / Rs
    if fm >= 0:
        mode = fm
    elif q == 1:
        mode = 0
    elif iL < 0.0:
        mode = 3
    elif iL > ith:
        mode = 1
    else:
        mode = 2
    vd = 0.0; idd = 0.0
    if mode == 0:
        gs = 1.0 / Ron + 1.0 / Rs
        if iL < 0.0:
            gs += 1.0 / Rb          # reverse current shared with the body diode
        vx = Vin - iL / gs
        isw = iL
    elif mode == 3:
        # body diode of the MOSFET (Vfd = 0, Rd) returns the current to Vin
        vx = Vin - iL / (1.0 / Rb + 1.0 / Rs)
        isw = iL
    elif mode == 1:
        vx = (Vin / Rs - Vf / rd - iL) / (1.0 / Rs + 1.0 / rd)
        idd = (-vx - Vf) / rd
        vd = -vx
        isw = (Vin - vx) / Rs
    else:
        vx = Vin - Rs * iL
        isw = iL
    diL = (vx - rL * iL - vo) / L
    if mode == 2:
        diL = 0.0          # quasi-static branch, handled outside
    pin = Vin * isw
    pu = vo * vo / R + vo * ip
    pls = rL * iL * iL + rC * ic * ic + (Vin - vx) * isw + vd * idd
    return diL, ic / C, pin, pu, pls, mode


@njit(cache=True)
def sim_sw(pl, p, tr, lv, ld_on, ld_off, ld_amp, iL0, vC0, T, h, fs, win):
    """Switched simulation. win: (3,2) steady windows for window energies.

    Returns Ts-sampled traces and a diagnostics vector:
      0 Ein 1 Eu 2 Eloss 3 Estored(T)-Estored(0) 4 residual (J)
      5 iL peak 6 vo max 7 DCM time fraction 8 saturation fraction
      9 reverse high-side conduction time 10 IAE (substep trapezoid)
      11..14 window residuals (%), 15..17 window efficiencies
    """
    Ts = p[32]
    Tsw = 1.0 / fs
    n = int(round(T / Ts))
    Vin, L, C, R, rL, rC, Ron, rd, Vf, Rs = pl[0], pl[1], pl[2], pl[3], pl[4], pl[5], pl[6], pl[7], pl[8], pl[9]
    s = np.zeros(NS)
    iL = iL0; vC = vC0
    r0 = _ref(0.0, tr, lv)
    ic, vo = _cap(iL, vC, 0.0, R, rC)
    ctrl_init(p, s, r0, vo, iL, vC)
    tv = np.zeros(n + 1); rv = np.zeros(n + 1); vv = np.zeros(n + 1)
    iv = np.zeros(n + 1); dv = np.zeros(n + 1); dpre = np.zeros(n + 1)
    fdcm = np.zeros(n + 1)
    Ein = 0.0; Eu = 0.0; El = 0.0
    Es0 = 0.5 * L * iL * iL + 0.5 * C * vC * vC
    # window accumulators: Ein, Eu, El, Es_start, Es_end
    wacc = np.zeros((win.shape[0], 5))
    for w in range(win.shape[0]):
        wacc[w, 3] = np.nan
    ipk = iL; vmax = vo
    t_dcm = 0.0; t_rev = 0.0; nsat = 0
    iae = 0.0
    eps = 1e-13
    for k in range(n):
        t = k * Ts
        r = _ref(t, tr, lv)
        ip = ld_amp if (ld_on <= t < ld_off) else 0.0
        ic, vo = _cap(iL, vC, ip, R, rC)
        dnow, dsat, dt_ = ctrl_step(p, s, r, vo, iL, ic)
        tv[k] = t; rv[k] = r; vv[k] = vo; iv[k] = iL; dv[k] = dnow; dpre[k] = dt_
        if dsat != dt_:
            nsat += 1
        for w in range(win.shape[0]):
            if abs(t - win[w, 0]) < 0.5 * Ts:
                wacc[w, 3] = 0.5 * L * iL * iL + 0.5 * C * vC * vC
                wacc[w, 0] = Ein; wacc[w, 1] = Eu; wacc[w, 2] = El
        ton = dnow * Tsw
        tdcm_k = t_dcm
        tl = t
        tstop = t + Ts
        while tl < tstop - eps:
            nper = np.floor((tl + eps) / Tsw)
            ph = tl - nper * Tsw
            if ph < ton - eps:
                q = 1
                tedge = nper * Tsw + ton
            else:
                q = 0
                tedge = (nper + 1.0) * Tsw
            if tedge <= tl + eps:
                tedge = tl + Tsw
            dt = min(tstop - tl, tedge - tl, h)
            if dt <= eps:
                dt = tstop - tl
            e0 = abs(r - vo)
            a1, b1, p1, u1, l1, m1 = _sw_f(iL, vC, q, ip, pl, -1)
            if m1 == 2:
                # quasi-static leakage branch: iL = (Vin - vo)/(Rs + rL)
                a2, b2, p2, u2, l2, m2 = _sw_f(iL, vC + 0.5 * dt * b1, q, ip, pl, m1)
                a3, b3, p3, u3, l3, m3 = _sw_f(iL, vC + 0.5 * dt * b2, q, ip, pl, m1)
                a4, b4, p4, u4, l4, m4 = _sw_f(iL, vC + dt * b3, q, ip, pl, m1)
                vC += dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
                Ein += dt / 6 * (p1 + 2 * p2 + 2 * p3 + p4)
                Eu += dt / 6 * (u1 + 2 * u2 + 2 * u3 + u4)
                El += dt / 6 * (l1 + 2 * l2 + 2 * l3 + l4)
                ic, vo = _cap(iL, vC, ip, R, rC)
                inew = (Vin - vo) / (Rs + rL)
                # stored-energy jump of the 1-ns branch is booked as loss
                El += 0.5 * L * (iL * iL - inew * inew)
                iL = inew
                t_dcm += dt
            else:
                a2, b2, p2, u2, l2, m2 = _sw_f(iL + 0.5 * dt * a1, vC + 0.5 * dt * b1, q, ip, pl, m1)
                a3, b3, p3, u3, l3, m3 = _sw_f(iL + 0.5 * dt * a2, vC + 0.5 * dt * b2, q, ip, pl, m1)
                a4, b4, p4, u4, l4, m4 = _sw_f(iL + dt * a3, vC + dt * b3, q, ip, pl, m1)
                iL += dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
                vC += dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
                Ein += dt / 6 * (p1 + 2 * p2 + 2 * p3 + p4)
                Eu += dt / 6 * (u1 + 2 * u2 + 2 * u3 + u4)
                El += dt / 6 * (l1 + 2 * l2 + 2 * l3 + l4)
                if q == 0 and ((m1 == 1 and iL < (Vin + Vf) / Rs) or (m1 == 3 and iL > 0.0)):
                    # diode blocks inside the step: jump to the leakage branch
                    ic, vo = _cap(iL, vC, ip, R, rC)
                    inew = (Vin - vo) / (Rs + rL)
                    El += 0.5 * L * (iL * iL - inew * inew)
                    iL = inew
                if iL < 0.0:
                    t_rev += dt
            ic, vo = _cap(iL, vC, ip, R, rC)
            iae += 0.5 * (e0 + abs(r - vo)) * dt
            if iL > ipk:
                ipk = iL
            if vo > vmax:
                vmax = vo
            tl += dt
        fdcm[k] = (t_dcm - tdcm_k) / Ts
        for w in range(win.shape[0]):
            if abs(tstop - win[w, 1]) < 0.5 * Ts:
                Es = 0.5 * L * iL * iL + 0.5 * C * vC * vC
                wacc[w, 0] = Ein - wacc[w, 0]
                wacc[w, 1] = Eu - wacc[w, 1]
                wacc[w, 2] = El - wacc[w, 2]
                wacc[w, 4] = Es - wacc[w, 3]
        if not (abs(vC) < 1e4 and abs(iL) < 1e4):
            break
    t = n * Ts
    r = _ref(t, tr, lv)
    ic, vo = _cap(iL, vC, 0.0, R, rC)
    tv[n] = t; rv[n] = r; vv[n] = vo; iv[n] = iL; dv[n] = s[4]; dpre[n] = s[5]
    dEs = 0.5 * L * iL * iL + 0.5 * C * vC * vC - Es0
    res = Ein - Eu - El - dEs
    nw = win.shape[0]
    diag = np.zeros(11 + 2 * (nw + 1))
    diag[0] = Ein; diag[1] = Eu; diag[2] = El; diag[3] = dEs; diag[4] = res
    diag[5] = ipk; diag[6] = vmax; diag[7] = t_dcm / T; diag[8] = nsat / n
    diag[9] = t_rev; diag[10] = iae
    diag[11] = 100 * res / Ein
    diag[11 + nw + 1] = (Eu + dEs * 0) / Ein
    for w in range(nw):
        rw = wacc[w, 0] - wacc[w, 1] - wacc[w, 2] - wacc[w, 4]
        diag[12 + w] = 100 * rw / wacc[w, 0]
        diag[13 + nw + w] = wacc[w, 1] / wacc[w, 0]
    return tv, rv, vv, iv, dv, dpre, diag, fdcm
