"""Matched-seed campaign: 6 global operators + refineBO, for every family.

usage: python scripts/cmo_campaign.py [family ...]
Writes results/cmo/campaign/<family>/<method>_<seed>.json (resumable).
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cmo.algos import METHODS, run_method, best_of, B_GLOBAL  # noqa: E402
from cmo.designs import FAMILIES, dec_cmo, XI0  # noqa: E402
from cmo.evaluate import Evaluator  # noqa: E402

SEEDS = [2026091200 + r for r in range(1, 31)]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'cmo', 'campaign')

_E = {}


def _job(args):
    fam, method, seed = args
    path = os.path.join(OUT, fam, f'{method}_{seed}.json')
    if os.path.exists(path):
        return path
    if fam not in _E:
        dec, lb, ub, x0 = FAMILIES[fam]()
        _E[fam] = (Evaluator(dec, dec_cmo(XI0)), lb, ub, x0)
    E, lb, ub, x0 = _E[fam]
    t = time.time()
    h = run_method(method, E, lb, ub, x0, seed)
    b, bg = best_of(h), best_of(h[:B_GLOBAL])
    curve, cur = [], None
    for r in h:
        if r.feasible and (cur is None or r.J < cur):
            cur = r.J
        curve.append(cur)
    rec = dict(family=fam, method=method, seed=seed, wall=time.time() - t,
               x=b.x.tolist(), J=b.J, g=b.g.tolist(), viol=b.viol, feasible=b.feasible,
               ripple=b.extra.get('ripple'), sat=b.extra.get('sat'), clip=b.extra.get('clip'),
               ratios=np.asarray(b.extra.get('ratios')).tolist(),
               xg=bg.x.tolist(), Jg=bg.J, feasible_g=bg.feasible, viol_g=bg.viol,
               n_feasible=int(sum(r.feasible for r in h)), curve=curve,
               stage=b.stage, n_eval=b.n_eval)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(rec, f)
    return path


if __name__ == '__main__':
    fams = sys.argv[1:] or list(FAMILIES)
    jobs = [(f, m, s) for f in fams for s in SEEDS for m in METHODS]
    with Pool(4) as pool:
        for i, p in enumerate(pool.imap_unordered(_job, jobs), 1):
            if i % 20 == 0:
                print(f'{i}/{len(jobs)} {p}', flush=True)
    print('done')
