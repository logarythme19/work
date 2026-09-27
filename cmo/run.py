"""Python-level wrappers: parameter vectors, scenarios, and replays."""
from __future__ import annotations

import numpy as np

from .params import PLANT, SAMP, MISSION_T, MISSION_STEPS, MISSION_LEVELS, MISSION_LOAD, MISSION_WINDOWS
from .ctrl import base_vector
from . import sim


def plant_vec(P=PLANT):
    return np.array([P.Vin, P.L, P.C, P.R, P.rL, P.rC, P.Ron, P.rd, P.Vf, P.Rs, P.Rbody])


def pidf_vector(theta, P=PLANT, S=SAMP):
    p = base_vector(P, S)
    p[0] = 3
    p[1:8] = theta[[0, 1, 2, 3, 4, 5, 6]]
    return p


def mission_sw(p, h=0.2e-6, P=PLANT, T=MISSION_T):
    tr = np.array(MISSION_STEPS); lv = np.array(MISSION_LEVELS)
    on, off, amp = MISSION_LOAD
    win = np.array(MISSION_WINDOWS)
    out = sim.sim_sw(plant_vec(P), p, tr, lv, on, off, amp, 2.4, 24.0, T, h, SAMP.fs, win)
    keys = ['t', 'r', 'vo', 'iL', 'd', 'dpre']
    res = dict(zip(keys, out[:6]))
    dg = out[6]
    res['fdcm'] = out[7]
    names = ['Ein', 'Eu', 'Eloss', 'dEs', 'res', 'ipk', 'vmax', 'dcm', 'sat', 'trev', 'IAE',
             'res_pct', 'res_w24', 'res_w1', 'res_w46', 'eta', 'eta_w24', 'eta_w1', 'eta_w46']
    res.update({k: float(v) for k, v in zip(names, dg)})
    return res
