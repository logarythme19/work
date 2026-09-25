"""Campaign statistics, freeze, switched validation and architecture comparison.

Outputs results/cmo/analysis.json and results/cmo/controllers.csv (the frozen
controllers, for matlab/run_cmo_validation.m).
"""
from __future__ import annotations

import glob
import json
import os
import sys
from itertools import combinations

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from cmo.algos import METHODS  # noqa: E402
from cmo.designs import FAMILIES, fixed_comparators, dec_cmo, XI0  # noqa: E402
from cmo.evaluate import Evaluator  # noqa: E402
from cmo.run import mission_sw  # noqa: E402
from cmo.params import MISSION_WINDOWS, MISSION_STEPS, IDS_RATING_A  # noqa: E402
from cmo.ripple import ripple_ratio  # noqa: E402
from cmo.lqi import lqi  # noqa: E402

OUT = os.path.join(ROOT, 'results', 'cmo')


def load():
    runs = {}
    for f in glob.glob(os.path.join(OUT, 'campaign', '*', '*.json')):
        r = json.load(open(f))
        runs.setdefault(r['family'], []).append(r)
    return runs


def score(r, key='J', feas='feasible', viol='viol'):
    """Deb-consistent scalar: J if feasible, 1e3 + violation otherwise."""
    return r[key] if r[feas] else 1e3 + r[viol]


def holm(p):
    k = sorted(p, key=p.get)
    m, run, out = len(k), 0.0, {}
    for i, key in enumerate(k):
        run = max(run, min(1.0, (m - i) * p[key]))
        out[key] = run
    return out


def family_stats(runs):
    by = {m: {r['seed']: r for r in runs if r['method'] == m} for m in METHODS}
    seeds = sorted(set.intersection(*[set(v) for v in by.values()]))
    S = np.array([[score(by[m][s]) for m in METHODS] for s in seeds])
    Sg = np.array([[score(by[m][s], 'Jg', 'feasible_g', 'viol_g') for m in METHODS] for s in seeds])
    res = {'n_seeds': len(seeds), 'methods': {}}
    for j, m in enumerate(METHODS):
        rr = [by[m][s] for s in seeds]
        Jf = np.array([r['J'] for r in rr if r['feasible']])
        res['methods'][m] = dict(
            feasible=int(len(Jf)), best=float(Jf.min()) if len(Jf) else None,
            median=float(np.median(Jf)) if len(Jf) else None,
            iqr=float(np.subtract(*np.percentile(Jf, [75, 25]))) if len(Jf) else None,
            worst=float(Jf.max()) if len(Jf) else None,
            feasible_global=int(sum(r['feasible_g'] for r in rr)),
            median_global=float(np.median([score(r, 'Jg', 'feasible_g', 'viol_g') for r in rr])),
            median_final=float(np.median(S[:, j])),
            viol_best=float(min(r['viol'] for r in rr)), viol_median=float(np.median([r['viol'] for r in rr])),
            ripple_med=float(np.median([r['ripple'] for r in rr])),
            ripple_max=float(np.max([r['ripple'] for r in rr])),
            refine_gain_pct=float(np.median(100 * (Sg[:, j] - S[:, j]) / Sg[:, j])),
            wall=float(np.mean([r['wall'] for r in rr])))
    if len(seeds) >= 5:
        fr = stats.friedmanchisquare(*S.T)
        ranks = np.mean(stats.rankdata(S, axis=1), axis=0)
        res['friedman'] = dict(Q=float(fr.statistic), p=float(fr.pvalue),
                               mean_ranks=dict(zip(METHODS, ranks.tolist())))
        raw = {}
        for a, b in combinations(range(len(METHODS)), 2):
            d = S[:, a] - S[:, b]
            raw[f'{METHODS[a]}|{METHODS[b]}'] = 1.0 if np.all(d == 0) else float(stats.wilcoxon(S[:, a], S[:, b]).pvalue)
        res['wilcoxon_holm'] = {k: dict(p=raw[k], p_holm=v) for k, v in holm(raw).items()}
    # convergence (median best-so-far feasible J)
    curves = []
    for m in METHODS:
        C = np.array([[np.nan if v is None else v for v in by[m][s]['curve']] for s in seeds])
        curves.append(np.nanmedian(np.where(np.isnan(C), np.inf, C), axis=0).tolist())
    res['curves'] = dict(zip(METHODS, curves))
    dk = lambda r: (0, r['J']) if r['feasible'] else (1, r['viol'])
    best = min(runs, key=dk)
    keys = ('method', 'seed', 'J', 'x', 'ripple', 'feasible', 'viol')
    res['best'] = {k: best[k] for k in keys}
    res['branch_best'] = {m: {k: b[k] for k in keys} for m in METHODS
                          for b in [min((r for r in runs if r['method'] == m), key=dk)]}
    return res, S


def sw_metrics(p):
    s = mission_sw(p)
    t = s['t']
    ess = []
    for a, b in MISSION_WINDOWS:
        m = (t >= a) & (t < b)
        ess.append(float(np.mean(s['r'][m] - s['vo'][m])))
    t1, t2 = MISSION_STEPS
    os46 = float(np.max(s['vo'][t >= t2]) - 46)
    us1 = float(np.max(-(s['vo'][(t >= t1) & (t < t2)] - 1.0)))
    m1 = (t >= MISSION_WINDOWS[1][0]) & (t < MISSION_WINDOWS[1][1])
    dpp = max(float(np.ptp(s['d'][(t >= a) & (t < b)])) for a, b in MISSION_WINDOWS)
    mload = (t >= 0.054) & (t < 0.108)
    ld = float(np.max(np.abs(s['r'][mload] - s['vo'][mload])))
    # mission specification ratios (<= 1 means compliant)
    spec = dict(overshoot=max(os46, us1, 0.0) / 1.5, load_dev=ld / 0.5, current=s['ipk'] / 23.0,
                steady=max(abs(e) for e in ess) / 0.025, duty_pp=dpp / 0.20)
    return dict(IAE=s['IAE'], ipk=s['ipk'], ipk_ratio=s['ipk'] / IDS_RATING_A, vmax=s['vmax'],
                os46=os46, ess24=ess[0], ess1=ess[1], ess46=ess[2],
                ess_worst=float(max(abs(e) for e in ess)), eta=s['eta'], LI=s['Eloss'] / s['Eu'],
                Eu=s['Eu'], Eloss=s['Eloss'], res_pct=s['res_pct'],
                res_windows=[s['res_w24'], s['res_w1'], s['res_w46']],
                eta_w1=s['eta_w1'], dcm_w1=float(np.mean(s['fdcm'][m1])), sat=s['sat'],
                trev_ms=1e3 * s['trev'], ripple=ripple_ratio(p), us1=us1, load_dev=ld, dpp=dpp,
                spec=spec, spec_util=max(spec.values()), spec_pass=int(sum(v <= 1 for v in spec.values())),
                spec_ok=bool(max(spec.values()) <= 1.0))


def main():
    runs = load()
    out = {'families': {}, 'switched': {}}
    fam_scores = {}
    for fam, rr in sorted(runs.items()):
        st, S = family_stats(rr)
        out['families'][fam] = st
        fam_scores[fam] = S
    # paired family ablations (same method, same seed): CMO vs the others
    abl = {}
    if 'CMO-LQI-PIDF' in fam_scores:
        A = fam_scores['CMO-LQI-PIDF'].ravel()
        for fam, S in fam_scores.items():
            if fam == 'CMO-LQI-PIDF' or S.shape != fam_scores['CMO-LQI-PIDF'].shape:
                continue
            B = S.ravel()
            d = A - B
            p = 1.0 if np.all(d == 0) else float(stats.wilcoxon(A, B).pvalue)
            abl[fam] = dict(n=len(A), cmo_better=int(np.sum(d < 0)), ties=int(np.sum(d == 0)),
                            other_better=int(np.sum(d > 0)), p=p,
                            median_cmo=float(np.median(A)), median_other=float(np.median(B)),
                            feasible_cmo=int(np.sum(A < 1e3)), feasible_other=int(np.sum(B < 1e3)))
        out['ablation_vs_CMO'] = abl
    # freeze: best feasible record of each family, plus best of each CMO branch
    ctrls = {}
    for fam, st in out['families'].items():
        dec = FAMILIES[fam]()[0]
        ctrls[f'{fam}*'] = dec(np.array(st['best']['x']))
    for m, b in out['families'].get('CMO-LQI-PIDF', {}).get('branch_best', {}).items():
        ctrls[f'CMO[{m}]'] = dec_cmo(np.array(b['x']))
    ctrls['Bryson anchor (xi0)'] = dec_cmo(XI0)
    ctrls.update(fixed_comparators())
    for k, p in ctrls.items():
        out['switched'][k] = sw_metrics(p)
        out['switched'][k]['p'] = p.tolist()
    # averaged-evaluator score of every frozen controller on the same scoreboard
    E = Evaluator(dec_cmo, dec_cmo(XI0))
    for k, p in ctrls.items():
        r = E.score(p)
        out['switched'][k]['avg_J'] = r.J
        out['switched'][k]['avg_viol'] = r.viol
        out['switched'][k]['avg_ratios'] = np.asarray(r.extra.get('ratios')).tolist() if 'ratios' in r.extra else None
    if 'CMO-LQI-PIDF' in out['families']:
        D = lqi(np.array(out['families']['CMO-LQI-PIDF']['best']['x']))
        out['retained_lqi'] = dict(K=D.K.tolist(), theta=D.theta.tolist(), poles=[[z.real, z.imag] for z in D.poles],
                                   care_res_rel=D.care_res_rel)
    with open(os.path.join(OUT, 'analysis.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    with open(os.path.join(OUT, 'controllers.csv'), 'w') as f:
        f.write('name,' + ','.join(f'p{i + 1}' for i in range(36)) + '\n')
        for k, p in ctrls.items():
            f.write('"' + k + '",' + ','.join(repr(float(v)) for v in p) + '\n')
    for k, v in out['switched'].items():
        print(f"{k:28s} IAE={v['IAE']:.5f} ipk={v['ipk']:.2f} os46={v['os46']:+.3f} ess={v['ess_worst']:.4f} "
              f"eta={100 * v['eta']:.4f} LI={v['LI']:.6f} sat={v['sat']:.3f} util={v['spec_util']:.2f} pass={v['spec_pass']}/5 avgJ={v['avg_J']:.4f} viol={v['avg_viol']:.3f}")


if __name__ == '__main__':
    main()
