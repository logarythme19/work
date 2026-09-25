"""Analytical switching-ripple bound (constraint g37).

The averaged evaluator contains no switching ripple. In the sampled loop the
controller runs at Ts = 1/(2 fs): the carrier-frequency output ripple is seen
exactly at the controller Nyquist frequency z = -1, where it appears as a
sample-to-sample alternation of amplitude at most dVo_pp/2 (CCM, triangular
inductor ripple, ESR + capacitive terms). The duty excursion it injects is

    rho_rip = max_{V in {1, 24, 46}} |C_fb(z = -1)| * dVo_pp(V) / 2,

with C_fb the discrete feedback path from -v_o to d. The bound rho_rip <= 0.20
keeps the injected alternation below one fifth of the carrier amplitude.
For state feedback the same construction uses the measured-state ripples.
"""
from __future__ import annotations

import numpy as np

from .params import PLANT, SAMP
from .plant import output_ripple_pp, equilibrium

LEVELS = (1.0, 24.0, 46.0)


def c_nyquist(p):
    code = int(p[0])
    Ts = p[32]
    if code == 1:
        return abs(p[1] - p[2] * Ts / 2)      # forward-Euler integrator Ts/(z-1) = -Ts/2 at z=-1
    if code == 2:
        return abs(p[1] - p[2] * Ts / 2 + 2 * p[3] / Ts)
    if code == 3:
        Kp, Kd, N = p[1], p[3], p[4]
        phi = np.exp(-N * Ts)
        return abs(Kp - p[2] * Ts / 2 + Kd * N * 2.0 / (1.0 + phi))
    return np.nan


def ripple_ratio(p, P=PLANT):
    code = int(p[0])
    out = 0.0
    for V in LEVELS:
        dvo = output_ripple_pp(V, P)
        if code in (1, 2, 3):
            val = c_nyquist(p) * dvo / 2
        elif code in (4, 6):
            _, d = equilibrium(V, 0.0, P)
            dI = V * (1 - d) / (P.L * SAMP.fs)
            val = abs(p[8]) * dI / 2 + abs(p[9]) * dvo / 2
        else:
            return np.nan
        out = max(out, val)
    return float(out)
