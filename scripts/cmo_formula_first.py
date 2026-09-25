"""Formula-first chain, Fig. 1 -> Fig. 2 -> exact LQI <-> 2-DOF PIDF equivalence.

Writes results/cmo/formula_first.json with every analytical number of the paper.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cmo.params import PLANT, SAMP, MISSION_STEPS, MISSION_LEVELS, MISSION_T  # noqa: E402
from cmo.plant import ccm_screen, linearize, Gvd_coeffs, Gvg_coeffs, output_ripple_pp  # noqa: E402
from cmo.lqi import lqi, XI0, realization_audit, return_difference, augmented  # noqa: E402
from cmo.lmi import certificate  # noqa: E402
from cmo.equiv import sim_cont  # noqa: E402
from cmo.run import plant_vec  # noqa: E402
from cmo.ripple import ripple_ratio  # noqa: E402
from cmo.designs import dec_cmo  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'cmo')


def jl(x):
    if isinstance(x, dict):
        return {k: jl(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jl(v) for v in x]
    if isinstance(x, np.ndarray):
        if np.iscomplexobj(x):
            return [[float(v.real), float(v.imag)] for v in x.ravel()]
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return float(x)
    if isinstance(x, complex):
        return [x.real, x.imag]
    return x


def equivalence_tests(D):
    pl = plant_vec()
    tr = np.array(MISSION_STEPS); lv = np.array(MISSION_LEVELS)
    K, th = D.K, D.theta

    def pair(amp=0.0, plv=pl, Rn=PLANT.R):
        a = sim_cont(0, plv, Rn, K, th, th[6], tr, lv, 0.054, 0.090, amp, 2.4, 24.0, MISSION_T, 0.05e-6, SAMP.dmin, SAMP.dmax, 20000)
        b = sim_cont(1, plv, Rn, K, th, th[6], tr, lv, 0.054, 0.090, amp, 2.4, 24.0, MISSION_T, 0.05e-6, SAMP.dmin, SAMP.dmax, 20000)
        sat = float(np.mean((a[:, 3] <= SAMP.dmin + 1e-12) | (a[:, 3] >= SAMP.dmax - 1e-12)))
        return dict(max_dvo=float(np.max(np.abs(a[:, 1] - b[:, 1]))),
                    max_diL=float(np.max(np.abs(a[:, 2] - b[:, 2]))),
                    max_dd=float(np.max(np.abs(a[:, 3] - b[:, 3]))),
                    sat_fraction=sat, iL_min=float(a[:, 2].min()))
    out = {'nominal_ip0': pair()}
    out['load_pulse_0p5A'] = pair(0.5)
    p2 = pl.copy(); p2[3] = 1.2 * PLANT.R
    out['R_plus20pct'] = pair(0.0, p2)
    p3 = pl.copy(); p3[1] = 0.8 * PLANT.L; p3[4] = 2 * PLANT.rL; p3[6] = 2 * PLANT.Ron
    out['L_minus20_rL_Ron_x2'] = pair(0.0, p3)
    p4 = pl.copy(); p4[2] = 1.2 * PLANT.C
    out['C_plus20pct'] = pair(0.0, p4)
    return out


def main():
    res = {}
    res['ccm'] = [ccm_screen(V) for V in (1.0, 24.0, 46.0)]
    tfs = []
    for V in (1.0, 24.0, 46.0):
        lin = linearize(V)
        n, d = Gvd_coeffs(V); ng, _ = Gvg_coeffs(V)
        tfs.append(dict(Vo=V, IL0=lin.IL0, d0=lin.d0, rho0=lin.rho0, Ed=lin.Ed,
                        den_monic=(d / d[0]).tolist(), Gvd_num=(n / d[0]).tolist(),
                        Gvg_num=(ng / d[0]).tolist(),
                        w0=float(np.sqrt(d[2] / d[0])), zeta0=float(d[1] / d[0] / (2 * np.sqrt(d[2] / d[0]))),
                        ripple_pp=output_ripple_pp(V)))
    res['transfer_functions'] = tfs
    lin = linearize(24.0)
    Aa, Ba, _ = augmented()
    ctrb = np.linalg.matrix_rank(np.hstack([lin.Bd, lin.A @ lin.Bd]))
    obsv = np.linalg.matrix_rank(np.vstack([lin.Cy, lin.Cy @ lin.A]))
    ctrb_a = np.linalg.matrix_rank(np.hstack([Ba, Aa @ Ba, Aa @ Aa @ Ba]))
    res['linear_24V'] = dict(A=lin.A, Bd=lin.Bd.ravel(), Cy=lin.Cy.ravel(), Dp=lin.Dp,
                             gamma=lin.gamma, rank_ctrb=ctrb, rank_obsv=obsv, rank_ctrb_aug=ctrb_a)
    D = lqi(XI0)
    coef, rel = realization_audit(D.K)
    rdmin, pm = return_difference(D)
    res['anchor'] = dict(xi0=XI0, Q=np.diag(D.Q), Rd=float(D.Rd[0, 0]), K=D.K,
                         care_res=D.care_res, care_res_rel=D.care_res_rel,
                         poles=D.poles, theta=D.theta, coef_residual=coef, freq_rel_error=rel,
                         min_return_difference=rdmin, phase_margin_deg=pm,
                         NC=1 / (PLANT.rC * PLANT.C), nyquist=np.pi / SAMP.Ts,
                         ripple_ratio=ripple_ratio(dec_cmo(XI0)))
    res['anchor_lmi'] = certificate(D.K)
    res['equivalence'] = equivalence_tests(D)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'formula_first.json'), 'w') as f:
        json.dump(jl(res), f, indent=1)
    print(json.dumps(jl({k: res[k] for k in ('anchor', 'equivalence')}), indent=1)[:4000])


if __name__ == '__main__':
    main()
