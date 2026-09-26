"""Noyau de simulation compile (numba) de la mission 24 -> 1 -> 46 V.

MODELE COMMUTE, pas moyenne. Ce choix n'est pas un raffinement optionnel:
l'espace de travail du modele Simulink de reference documente explicitement que

    "The 1-V setpoint sits 3 % from the CCM/DCM boundary at the nominal 10 ohm
     load, and the 24 -> 1 V step forces about 3.2 ms of discontinuous
     operation."

Un evaluateur moyenne a conduction continue est donc structurellement faux
exactement la ou vit la contrainte g4 (|e1| <= 0.03 V au palier 1 V). La diode
de roue libre est ici modelisee avec son blocage reel, ce qui produit la
conduction discontinue sans modele DCM separe.

Trois etats de conduction:
    1. interrupteur ON            : L diL/dt = Vin - vo - (Ron + rL) iL
    2. OFF, diode passante        : L diL/dt = -Vf - (rd + rL) iL - vo
    3. OFF, diode bloquee (iL=0)  : diL/dt = 0                        (DCM)

Le rapport cyclique est echantillonne a Ts avec un retard d'un echantillon,
tandis que la porteuse PWM est continue a 50 kHz. Les fronts de commutation
sont traites par DECOUPAGE EXACT du pas d'integration: sans cela, un pas de
1 us sur une periode de 20 us ne donnerait que 5 % de resolution en rapport
cyclique, ce qui rendrait les bornes [0.02, 0.98] illusoires.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, inline='always')
def _cap(iL, vC, ip, R, rC, C):
    """Courant condensateur, tension de sortie et derivee de vC."""
    ic = (iL - vC / R - ip) / (1.0 + rC / R)
    vo = vC + rC * ic
    return ic, vo, ic / C


@njit(cache=True, inline='always')
def _diL(iL, vo, mode, Vin, L, rL, Ron, rd, Vf):
    if mode == 0:        # interrupteur ferme
        return (Vin - vo - (Ron + rL) * iL) / L
    elif mode == 1:      # diode de roue libre passante
        return (-Vf - (rd + rL) * iL - vo) / L
    else:                # diode bloquee: courant maintenu a zero
        return 0.0


@njit(cache=True)
def simulate(
    Vin, L, C, R, rL, rC, Ron, rd, Vf,
    v0, v1, v2, t1, t2, t_end,
    ip_amp, ip_t0, ip_t1, iL0, vC0,
    Ts, h, Tsw, duty_min, duty_max, delay_n,
    Vn, Dn,
    Rn, rLn, rdn, Ronn, Vinn,
    kp, ki, kd, wf, bsp, csp, Kb,
    zi, pi_, phii, betai, gaini,
    zd, pd, phid, betad, gaind,
    w24a, w24b, w1a, w1b, w46a, w46b,
):
    n_ctrl = int(round(t_end / Ts))

    iL = iL0
    vC = vC0

    ni = zi.shape[0]
    nd = zd.shape[0]
    wi = np.zeros(ni)
    wd = np.zeros(nd)

    xI = 0.0
    xF = 0.0
    phiF = np.exp(-wf * Ts)

    dff0 = (Rn + rLn + rdn) * v0 / (Rn * Vinn - (Ronn - rdn) * v0)
    # Registre a decalage de longueur delay_n: d applique a l'instant k est
    # celui calcule a l'instant k - delay_n.
    dbuf = np.full(max(delay_n, 1), dff0)

    int_err = 0.0
    int_iL2 = 0.0
    sum_dd = 0.0
    n_duty = 0
    d_prev = dff0

    Imax = 0.0
    Vmax = -1e30
    Vmin = 1e30

    s24 = 0.0; n24 = 0
    s1 = 0.0;  n1 = 0
    s46 = 0.0; n46 = 0

    n_dcm = 0
    n_sub_tot = 0
    diverged = 0

    prev_err = 0.0
    prev_iL2 = 0.0
    have_prev = 0

    diode_on = 0

    for k in range(n_ctrl):
        t = k * Ts

        if t < t1:
            r = v0
        elif t < t2:
            r = v1
        else:
            r = v2

        ip = ip_amp if (ip_t0 <= t < ip_t1) else 0.0

        _, vo, _ = _cap(iL, vC, ip, R, rC, C)

        if not (vo == vo) or abs(vo) > 1e5 or abs(iL) > 1e4:
            diverged = 1
            break

        # ----------------------------------------------- correcteur FO-PIDF
        en = (r - vo) / Vn
        eb = (bsp * r - vo) / Vn
        ec = (csp * r - vo) / Vn

        y = en
        for i in range(ni):
            yi = y + (zi[i] - pi_[i]) * wi[i]
            wi[i] = phii[i] * wi[i] + betai[i] * y
            y = yi
        q = gaini * y

        xF = phiF * xF + (1.0 - phiF) * ec
        y = xF
        for i in range(nd):
            yi = y + (zd[i] - pd[i]) * wd[i]
            wd[i] = phid[i] * wd[i] + betad[i] * y
            y = yi
        dterm = gaind * y

        u = kp * eb + xI + kd * dterm

        # Feedforward NOMINAL fige (par.1), jamais les parametres du plant simule.
        dff = (Rn + rLn + rdn) * r / (Rn * Vinn - (Ronn - rdn) * r)
        d_tilde = dff + Dn * u
        d_sat = min(max(d_tilde, duty_min), duty_max)

        xI = xI + Ts * (ki * q + Kb * (d_sat - d_tilde) / Dn)

        d_app = dbuf[0]
        for _s in range(delay_n - 1):
            dbuf[_s] = dbuf[_s + 1]
        dbuf[delay_n - 1] = d_sat

        sum_dd += abs(d_sat - d_prev)
        d_prev = d_sat
        n_duty += 1

        # ------------------------------------------ quadrature grille 10 us
        denom = r if r > 1.0 else 1.0
        cur_err = abs(r - vo) / denom
        cur_iL2 = iL * iL
        if have_prev == 1:
            int_err += 0.5 * (prev_err + cur_err) * Ts
            int_iL2 += 0.5 * (prev_iL2 + cur_iL2) * Ts
        prev_err = cur_err
        prev_iL2 = cur_iL2
        have_prev = 1

        if w24a <= t <= w24b:
            s24 += (r - vo); n24 += 1
        if w1a <= t <= w1b:
            s1 += (r - vo);  n1 += 1
        if w46a <= t <= w46b:
            s46 += (r - vo); n46 += 1

        # ------------------------------ integration commutee sur [t, t+Ts]
        t_local = t
        t_stop = t + Ts
        d_on_time = d_app * Tsw

        eps = 1e-13
        while t_local < t_stop - eps:
            n_per = np.floor((t_local + eps) / Tsw)
            phase = t_local - n_per * Tsw
            if phase < 0.0:
                phase = 0.0

            if phase < d_on_time - eps:
                sw_on = 1
                t_edge = n_per * Tsw + d_on_time
            else:
                sw_on = 0
                t_edge = (n_per + 1.0) * Tsw

            # Garde-fou: un front deja atteint ne doit jamais produire un pas
            # nul, sinon la periode de commande serait amputee.
            if t_edge <= t_local + eps:
                t_edge = t_local + Tsw

            dt = t_stop - t_local
            if t_edge - t_local < dt:
                dt = t_edge - t_local
            if h < dt:
                dt = h
            if dt <= eps:
                dt = t_stop - t_local
                if dt <= 0.0:
                    break

            # etat de conduction
            if sw_on == 1:
                mode = 0
                diode_on = 0
            else:
                _, vo_r, _ = _cap(iL, vC, ip, R, rC, C)
                if iL > 0.0:
                    mode = 1
                    diode_on = 1
                else:
                    # la diode ne conduit que si le courant tente de croitre
                    if (-Vf - vo_r) / L > 0.0:
                        mode = 1
                        diode_on = 1
                    else:
                        mode = 2
                        diode_on = 0
                        iL = 0.0
                        n_dcm += 1

            ic0, vo0, _ = _cap(iL, vC, ip, R, rC, C)
            if vo0 > Vmax:
                Vmax = vo0
            if vo0 < Vmin:
                Vmin = vo0
            if iL > Imax:
                Imax = iL
            n_sub_tot += 1

            k1i = _diL(iL, vo0, mode, Vin, L, rL, Ron, rd, Vf)
            _, _, k1v = _cap(iL, vC, ip, R, rC, C)

            i2 = iL + 0.5 * dt * k1i
            v2_ = vC + 0.5 * dt * k1v
            _, vo2, k2v = _cap(i2, v2_, ip, R, rC, C)
            k2i = _diL(i2, vo2, mode, Vin, L, rL, Ron, rd, Vf)

            i3 = iL + 0.5 * dt * k2i
            v3 = vC + 0.5 * dt * k2v
            _, vo3, k3v = _cap(i3, v3, ip, R, rC, C)
            k3i = _diL(i3, vo3, mode, Vin, L, rL, Ron, rd, Vf)

            i4 = iL + dt * k3i
            v4 = vC + dt * k3v
            _, vo4, k4v = _cap(i4, v4, ip, R, rC, C)
            k4i = _diL(i4, vo4, mode, Vin, L, rL, Ron, rd, Vf)

            iL = iL + (dt / 6.0) * (k1i + 2.0 * k2i + 2.0 * k3i + k4i)
            vC = vC + (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)

            # blocage de la diode: le courant d'inductance ne peut pas
            # devenir negatif tant que l'interrupteur est ouvert
            if mode == 1 and iL < 0.0:
                iL = 0.0
            if mode == 2:
                iL = 0.0

            t_local += dt

        if not (iL == iL) or not (vC == vC):
            diverged = 1
            break

    if diverged == 0:
        _, vo_r, _ = _cap(iL, vC, 0.0, R, rC, C)
        if vo_r > Vmax:
            Vmax = vo_r
        if vo_r < Vmin:
            Vmin = vo_r

    e24 = s24 / n24 if n24 > 0 else 1e30
    e1 = s1 / n1 if n1 > 0 else 1e30
    e46 = s46 / n46 if n46 > 0 else 1e30

    return (int_err, int_iL2, sum_dd, n_duty,
            Imax, Vmax, Vmin, e24, e1, e46,
            n_dcm, n_sub_tot, diverged)
