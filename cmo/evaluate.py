"""Averaged six-scenario evaluator, objective J_s and the 37 inequalities.

Scenarios (Table 5): S0 composite 24->1->46 V (steps 30/60 ms) with a 0.5 A
pulse 15-25 ms; S1a 24->1; S1b 1->46; S1c 46->24; S2a / S2b 0.5 A pulse at
24 / 46 V. Initial states are the analytical equilibria before the first event.

Per scenario: r_s = [rho_loss, rho_IAE, rho_ISE, rho_ITAEg, rho_ITAEev,
rho_RMSE, rho_MAE] (ratios to the normalization record xi0),
J_s = mean over the 42 entries (Eq. (42)).
Constraints (absolute engineering specifications, fixed before any search), g <= 0:
  g1 overshoot beyond the new level, in the step direction  <= 1.5 V
  g2 largest |e| from the load-pulse onset to the next event <= 0.5 V
  g3 peak inductor current <= 23 A (IRF540N continuous rating at 100 C)
  g4 largest final-window |mean error| <= 25 mV
  g5 |Eu/Eu0 - 1| <= 0.5 %   (same useful task as the record xi0)
  g6 loss intensity rho_loss <= 1.005
  g37 switched screen: duty peak-to-peak, last 4 ms of each level <= 0.20
  g38 switched screen: |mean error|, last 4 ms of each level       <= 25 mV
  g39 switched screen: peak inductor current                       <= 23 A
The switched screen (sw_screen) runs the switched model of cmo/sim.py over the
transition sequence 24 -> 1 -> 46 -> 24 V (34 ms, 1 us RK4). It replaces the analytical ripple
bound of the v6 protocol, which did not prevent ripple-induced limit cycles
(cmo/ripple.py keeps that bound as a reported diagnostic). The averaged peak
current in g3 is augmented by the worst half-ripple dI_pp/2 = 1.17 A (24 V).
A family that does not apply to a scenario (no step / no load pulse) is -1.
Feasible: sum(max(g, 0)) <= 1e-4 (declared acceptance band).
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

from .params import PLANT, SAMP
from .plant import equilibrium
from .run import plant_vec
from . import sim
from .ripple import ripple_ratio
from .plant import output_ripple_pp

RIPPLE_BOUND = 0.20
FEAS_TOL = 1e-4
ESS_TOL = 0.025
EU_BAND = 0.005
ESS_WIN = 2e-3
OS_LIM = 1.5
LOAD_LIM = 0.5
I_LIM = 23.0
HALF_RIPPLE = 0.5 * max(output_ripple_pp(V) / np.hypot(PLANT.rC, 1 / (8 * PLANT.C * SAMP.fs)) for V in (1.0, 24.0, 46.0))
SW_T = 10e-3
SW_H = 1e-6
DPP_LIM = 0.20
LI_TOL = 0.005


@dataclass(frozen=True)
class Scenario:
    name: str
    T: float
    tr: tuple
    lv: tuple
    load: tuple      # (t_on, t_off, amp)
    resets: tuple


SCENARIOS = (
    Scenario('S0', 0.090, (0.030, 0.060), (24.0, 1.0, 46.0), (0.015, 0.025, 0.5), (0.0, 0.015, 0.025, 0.030, 0.060)),
    Scenario('S1a', 0.040, (0.010,), (24.0, 1.0), (0.0, 0.0, 0.0), (0.0, 0.010)),
    Scenario('S1b', 0.040, (0.010,), (1.0, 46.0), (0.0, 0.0, 0.0), (0.0, 0.010)),
    Scenario('S1c', 0.040, (0.010,), (46.0, 24.0), (0.0, 0.0, 0.0), (0.0, 0.010)),
    Scenario('S2a', 0.040, (), (24.0,), (0.010, 0.020, 0.5), (0.0, 0.010, 0.020)),
    Scenario('S2b', 0.040, (), (46.0,), (0.010, 0.020, 0.5), (0.0, 0.010, 0.020)),
)


SCREEN_TR = (0.002, 0.014, 0.024)
SCREEN_LV = (24.0, 1.0, 46.0, 24.0)
SCREEN_T = 0.034
SCREEN_WIN = ((0.010, 0.014), (0.020, 0.024), (0.030, 0.034))


def sw_screen(p, P=PLANT):
    """Switched transition screen: 24 -> 1 -> 46 -> 24 V in 34 ms (1 us RK4).

    Returns (max duty p-p, max |mean error|, peak iL) over the last 4 ms of the
    1, 46 and 24 V segments. The large transitions matter: the loop can hold
    coexisting switching limit cycles, and the one reached after a large step
    is not the one reached from equilibrium.
    """
    IL, _ = equilibrium(24.0, 0.0, P)
    out = sim.sim_sw(plant_vec(P), p, np.array(SCREEN_TR), np.array(SCREEN_LV), 0.0, 0.0, 0.0,
                     IL, 24.0, SCREEN_T, SW_H, SAMP.fs, np.zeros((0, 2)))
    t, r, vo, d = out[0], out[1], out[2], out[4]
    if not np.all(np.isfinite(vo)):
        return np.inf, np.inf, np.inf
    dpp, ess = 0.0, 0.0
    for a, b in SCREEN_WIN:
        m = (t >= a) & (t < b)
        dpp = max(dpp, float(np.ptp(d[m])))
        ess = max(ess, abs(float(np.mean(r[m] - vo[m]))))
    return dpp, ess, float(out[6][5])


def scenario_metrics(p, sc: Scenario, nsub=5, P=PLANT):
    IL0, _ = equilibrium(sc.lv[0], 0.0, P)
    tv, rv, vv, iv, dv, dp, dg = sim.sim_avg(
        plant_vec(P), p, np.array(sc.tr, float), np.array(sc.lv, float),
        sc.load[0], sc.load[1], sc.load[2], IL0, sc.lv[0], sc.T, nsub)
    if not np.isfinite(dg[0]):
        return None
    e = rv - vv
    ae = np.abs(e)
    T = sc.T
    IAE = np.trapezoid(ae, tv)
    ISE = np.trapezoid(e * e, tv)
    ITAEg = np.trapezoid(tv * ae, tv)
    rs = np.array(sc.resets + (T,))
    ITAEev = 0.0
    ess = 0.0
    for a, b in zip(rs[:-1], rs[1:]):
        m = (tv >= a) & (tv <= b)
        ITAEev += np.trapezoid((tv[m] - a) * ae[m], tv[m])
        w = (tv >= b - ESS_WIN) & (tv < b)
        ess = max(ess, abs(np.mean(e[w])))
    os_ = -np.inf
    for k, ts in enumerate(sc.tr):
        sgn = np.sign(sc.lv[k + 1] - sc.lv[k])
        nxt = sc.tr[k + 1] if k + 1 < len(sc.tr) else T
        m = (tv >= ts) & (tv < nxt)
        os_ = max(os_, np.max(sgn * (vv[m] - sc.lv[k + 1])))
    ld = -np.inf
    if sc.load[2] > 0:
        nxt = min([b for b in rs if b > sc.load[1] + 1e-12] + [T])
        m = (tv >= sc.load[0]) & (tv < nxt)
        ld = np.max(ae[m])
    Ein, Eu, El = dg[0], dg[1], dg[2]
    return dict(IAE=IAE, ISE=ISE, ITAEg=ITAEg, ITAEev=ITAEev,
                RMSE=np.sqrt(ISE / T), MAE=IAE / T, ILpk=dg[5], ess=ess,
                Eu=Eu, El=El, LI=El / Eu, clip=dg[3], sat=dg[4], OS=os_, LD=ld)


@dataclass
class Record:
    x: np.ndarray
    J: float
    g: np.ndarray
    viol: float
    feasible: bool
    extra: dict = field(default_factory=dict)
    n_eval: int = 0
    stage: str = ''


class Evaluator:
    """Maps a decision vector to a Deb record, normalized by the record xi0."""

    CRIT = ('LI', 'IAE', 'ISE', 'ITAEg', 'ITAEev', 'RMSE', 'MAE')

    def __init__(self, decode, p_ref):
        self.decode = decode        # x -> controller parameter vector p
        self.ref = [scenario_metrics(p_ref, sc) for sc in SCENARIOS]

    def metrics(self, p):
        return [scenario_metrics(p, sc) for sc in SCENARIOS]

    def score(self, p, x=None):
        rip = ripple_ratio(p)
        ms = self.metrics(p)
        if any(m is None for m in ms):
            g = np.full(39, 1e3)
            dpp = ess_sw = ipk_sw = np.inf
            return Record(np.asarray(x), 1e3, g, float(np.sum(g)), False,
                          dict(ripple=rip, dpp=dpp, ess_sw=ess_sw, ipk_sw=ipk_sw))
        r = np.array([[m[c] / m0[c] for c in self.CRIT] for m, m0 in zip(ms, self.ref)])
        J = float(r.mean())
        g = []
        for m, m0 in zip(ms, self.ref):
            g += [m['OS'] / OS_LIM - 1 if np.isfinite(m['OS']) else -1.0,
                  m['LD'] / LOAD_LIM - 1 if np.isfinite(m['LD']) else -1.0,
                  (m['ILpk'] + HALF_RIPPLE) / I_LIM - 1, m['ess'] / ESS_TOL - 1,
                  abs(m['Eu'] / m0['Eu'] - 1) / EU_BAND - 1,
                  (m['LI'] / m0['LI'] - 1) / LI_TOL - 1]
        dpp, ess_sw, ipk_sw = sw_screen(p)
        g.append(dpp / DPP_LIM - 1)
        g.append(ess_sw / ESS_TOL - 1)
        g.append(ipk_sw / I_LIM - 1)
        g = np.array(g)
        viol = float(np.sum(np.maximum(g, 0.0)))
        return Record(np.asarray(x) if x is not None else None, J, g, viol,
                      viol <= FEAS_TOL,
                      dict(ripple=rip, dpp=dpp, ess_sw=ess_sw, ipk_sw=ipk_sw, ratios=r, clip=max(m['clip'] for m in ms),
                           sat=max(m['sat'] for m in ms)))

    def __call__(self, x):
        return self.score(self.decode(np.asarray(x, float)), x)
