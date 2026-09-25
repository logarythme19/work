"""Common sampled controller layer (numba), shared by both simulators.

Mirrors the "Common sampled controller" S-function of the reference .slx:
inputs Vref, V_o, I_L, I_c, V_in; output d_sat with ONE-SAMPLE delay; the
non-ideal nominal feedforward d_ff(r) of Eq. (7); total duty saturation; and
back-calculation anti-windup. Only the control law differs between codes.

Parameter vector p (float64, length NP):
  p[0]  code          1 PI, 2 PID, 3 PIDF-2DOF (also CMO-LQI-PIDF), 4 LQR,
                      5 LQG, 6 LQI (full state, measured iL)
  p[1:8]              law coefficients (see below)
  p[8:24]             matrices for LQR/LQG/LQI
  p[24:33]            nominal plant for feedforward: R rL rd Ron Vin Vf rC C Ts
  p[33:35]            dmin dmax
Law coefficients:
  PI   : Kp Ki Kb
  PID  : Kp Ki Kd b Kb           (backward-difference derivative on c*r - vo, c = 0)
  PIDF : Kp Ki Kd N b c Kb       (ZOH-exact first-order derivative filter)
  LQR  : p[8:10] = K1 K2
  LQG  : p[8:10] = K1 K2, p[10:12] = Lo, p[12:16] = Ad, p[16:18] = Bd, p[18:20] = Cd
  LQI  : p[8:11] = K1 K2 K3 (discrete), p[1] = Kb
State vector s (length NS): s[0] integral, s[1] filter / previous output,
  s[2:4] observer state, s[4] applied duty (delay register), s[5] d_pre.
"""
from __future__ import annotations

import numpy as np
from numba import njit

NP = 36
NS = 8


@njit(cache=True)
def feedforward(r, p):
    R = p[24]; rL = p[25]; rd = p[26]; Ron = p[27]; Vin = p[28]; Vf = p[29]
    IL = r / R
    return (r + (rL + rd) * IL + Vf) / (Vin - (Ron - rd) * IL + Vf)


@njit(cache=True)
def ctrl_init(p, s, r0, vo0, iL0, vC0):
    for i in range(NS):
        s[i] = 0.0
    code = int(p[0])
    d0 = feedforward(r0, p)
    s[4] = d0
    s[5] = d0
    if code == 3:
        # consistent initialization of the derivative filter:
        # x_F(0) = c r - v_C(0)  (Proposition 3: x_F tracks c r - v_C)
        s[1] = p[6] * r0 - vC0
    elif code == 2:
        s[1] = -vo0
    elif code == 5:
        s[2] = iL0
        s[3] = vC0


@njit(cache=True)
def ctrl_step(p, s, r, vo, iL, iC):
    """One sample. Returns (d_applied_now, d_sat_new, d_pre_new).

    d_applied_now is the value computed at the previous sample (one-sample
    delay of the S-function); the new value is stored for the next sample.
    """
    code = int(p[0])
    Ts = p[32]
    rC = p[30]
    R = p[24]
    dmin = p[33]; dmax = p[34]
    dff = feedforward(r, p)
    e = r - vo
    d_now = s[4]
    u = 0.0
    if code == 1:                      # PI
        Kp = p[1]; Ki = p[2]
        u = Kp * e + s[0]
    elif code == 2:                    # PID, backward difference on -vo
        Kp = p[1]; Kd = p[3]; b = p[4]
        yd = -vo
        u = Kp * (b * r - vo) + s[0] + Kd * (yd - s[1]) / Ts
        s[1] = yd
    elif code == 3:                    # PIDF 2-DOF
        Kp = p[1]; Kd = p[3]; N = p[4]; b = p[5]; c = p[6]
        u = Kp * (b * r - vo) + s[0] + Kd * N * (c * r - vo - s[1])
        phi = np.exp(-N * Ts)
        s[1] = phi * s[1] + (1.0 - phi) * (c * r - vo)
    elif code == 4 or code == 6:       # LQR / LQI, measured state
        vC = vo - rC * iC
        u = -p[8] * (iL - r / R) - p[9] * (vC - r)
        if code == 6:
            u += -p[10] * s[0]
    elif code == 5:                    # LQG, Kalman predictor on vo only
        xh0 = s[2]; xh1 = s[3]
        xe0 = r / R; xe1 = r
        yh = r + p[18] * (xh0 - xe0) + p[19] * (xh1 - xe1)
        nu = vo - yh
        u = -p[8] * (xh0 - xe0) - p[9] * (xh1 - xe1)
        # predictor update uses the duty that is actually applied now
        a0 = xh0 + p[10] * nu - xe0
        a1 = xh1 + p[11] * nu - xe1
        s[2] = xe0 + p[12] * a0 + p[13] * a1 + p[16] * (d_now - dff)
        s[3] = xe1 + p[14] * a0 + p[15] * a1 + p[17] * (d_now - dff)
    dt = dff + u
    ds = min(max(dt, dmin), dmax)
    # integral / anti-windup update (forward Euler, as in the S-function)
    if code == 1:
        s[0] += Ts * (p[2] * e + p[3] * (ds - dt))
    elif code == 2:
        s[0] += Ts * (p[2] * e + p[5] * (ds - dt))
    elif code == 3:
        s[0] += Ts * (p[2] * e + p[7] * (ds - dt))
    elif code == 6:
        K3 = abs(p[10])
        s[0] += Ts * (e + p[1] / K3 * (ds - dt))
    s[4] = ds
    s[5] = dt
    return d_now, ds, dt


def base_vector(P, S):
    p = np.zeros(NP)
    p[24:33] = [P.R, P.rL, P.rd, P.Ron, P.Vin, P.Vf, P.rC, P.C, S.Ts]
    p[33:35] = [S.dmin, S.dmax]
    return p
