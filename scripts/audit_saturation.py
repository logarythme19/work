#!/usr/bin/env python3
"""Livrable 2: audit de saturation des bornes d'ancrage, par composante.

Balaye la boite de decision v5 et mesure, pour chaque composante de theta, la
fraction de projections qui touchent chaque borne. Aucune simulation n'est
necessaire: la projection est algebrique, donc l'audit est exhaustif et peu
couteux. C'est precisement pour cela qu'il aurait du etre fait avant de
depenser 1800 evaluations.
"""

import argparse, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.stats import qmc

from isafo.folqi import synthesize, V5_LB, V5_UB, V5_NAMES, THETA_NAMES
from isafo.anchor import SaturationAudit, project, ANCHOR_LB, ANCHOR_UB


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-samples", type=int, default=4096)
    ap.add_argument("--seed", type=int, default=2026092401)
    ap.add_argument("--vref", type=float, default=24.0)
    ap.add_argument("-o", "--out", default="results/audit_saturation.json")
    args = ap.parse_args()

    sob = qmc.Sobol(d=7, scramble=True, seed=args.seed)
    m = int(np.ceil(np.log2(max(args.n_samples, 2))))
    X = qmc.scale(sob.random_base2(m), V5_LB, V5_UB)

    audit = SaturationAudit()
    raws = []
    for x in X:
        syn = synthesize(x, vref=args.vref)
        raws.append(syn.theta_raw)
        project(syn, audit=audit)

    raws = np.asarray(raws, dtype=float)
    rep = audit.to_dict()
    rep["n_echantillons"] = int(X.shape[0])
    rep["vref_synthese"] = args.vref
    rep["graine_sobol"] = args.seed

    # Distribution des composantes reellement produites par la synthese.
    rep["classification"] = audit.classify(raws)
    rep["distribution_theta_brut"] = {}
    for j, name in enumerate(THETA_NAMES):
        col = raws[:, j]
        fin = col[np.isfinite(col)]
        rep["distribution_theta_brut"][name] = {
            "n_fini": int(fin.size),
            "min": float(np.min(fin)) if fin.size else None,
            "p05": float(np.percentile(fin, 5)) if fin.size else None,
            "median": float(np.median(fin)) if fin.size else None,
            "p95": float(np.percentile(fin, 95)) if fin.size else None,
            "max": float(np.max(fin)) if fin.size else None,
        }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(rep, f, indent=2)

    print(f"Echantillons Sobol : {X.shape[0]}  (vref = {args.vref} V)")
    print(f"Projections rejetees (gain non positif) : {audit.n_rejected}"
          f"  ({audit.n_rejected/X.shape[0]*100:.2f} %)\n")
    print(f"{'composante':<12} {'lb':>8} {'ub':>8} {'%lb':>8} {'%ub':>8}  "
          f"{'p05':>9} {'mediane':>9} {'p95':>9}")
    lo, hi = audit.fractions()
    for j, name in enumerate(THETA_NAMES):
        d = rep["distribution_theta_brut"][name]
        f = lambda v: f"{v:9.3f}" if v is not None else "        -"
        flag = "  <-- ALARME" if (lo[j] > 0.20 or hi[j] > 0.20) else ""
        print(f"{name:<12} {ANCHOR_LB[j]:8.3f} {ANCHOR_UB[j]:8.3f} "
              f"{lo[j]*100:7.1f}% {hi[j]*100:7.1f}%  "
              f"{f(d['p05'])} {f(d['median'])} {f(d['p95'])}{flag}")

    print("\nNature des saturations (> 20 %, seuil du par.2) :")
    any_path = False
    for name, c in rep["classification"].items():
        if c["nature"] == "aucune":
            continue
        any_path |= (c["nature"] == "pathologique")
        print(f"  - {name:<12} {c['fraction_saturee']*100:5.1f} %  "
              f"etendue={c['etendue_theta_brut']:.3e}  -> {c['nature'].upper()}")
    if not any_path:
        print("\n  Aucune saturation PATHOLOGIQUE: toute borne atteinte l'est par une")
        print("  composante constante de la synthese, donc sans perte de liberte.")
    print(f"\n-> {args.out}")


if __name__ == "__main__":
    main()
