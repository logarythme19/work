"""Plant, sampling and mission constants decoded from the reference .slx.

Source: simulink/modelWorkspace.mxarray of
buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx (decoded with scipy.io).
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Plant:
    Vin: float = 48.0     # V, constant source (VIN_POLICY)
    L: float = 100e-6     # H
    C: float = 100e-6     # F
    R: float = 10.0       # ohm
    rL: float = 0.10      # ohm
    rC: float = 0.05      # ohm
    Ron: float = 0.05     # ohm, IRF540N static surrogate
    rd: float = 0.020     # ohm, STPST10H100SF affine law
    Vf: float = 0.42      # V
    Rs: float = 1e5       # ohm, resistive snubber across the switch (Cs = inf)
    Rbody: float = 0.01   # ohm, MOSFET body diode (Rd = 0.01, Vfd = 0 in the .slx)


@dataclass(frozen=True)
class Sampling:
    Ts: float = 10e-6     # s, controller sample time (Te)
    fs: float = 50e3      # Hz, PWM carrier
    dmin: float = 0.02
    dmax: float = 0.98


PLANT = Plant()
SAMP = Sampling()

# Validation mission of the .slx (Vref_ts, Iload_ts, Tsim)
MISSION_T = 0.324
MISSION_STEPS = (0.108, 0.216)
MISSION_LEVELS = (24.0, 1.0, 46.0)
MISSION_LOAD = (0.054, 0.090, 0.5)          # t_on, t_off, amplitude (A)
# steady windows: last 10 % of each level (used for steady error and energy)
MISSION_WINDOWS = ((0.0972, 0.108), (0.2052, 0.216), (0.3132, 0.324))
IDS_RATING_A = 33.0                          # IRF540N continuous drain, 25 C


def mission_ref(t):
    t1, t2 = MISSION_STEPS
    return np.where(t < t1, 24.0, np.where(t < t2, 1.0, 46.0))
