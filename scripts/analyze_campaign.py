#!/usr/bin/env python3
"""Livrables 3 et 4: resultats par graine et par methode, echecs compris, puis
analyse statistique conforme au par.6 (unite = la graine)."""

import argparse, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.algos import METHODS
from isafo.stats import compare_methods, wilson


def load(d):
    rows = []
    for p in sorted(glob.glob(os.path.join(d, "*.json"))):
        with open(p) as f:
            rows.append(json.load(f))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    rows = load(args.dir)
    if not rows:
        sys.exit(f"aucun resultat dans {args.dir}")
    seeds = sorted({r["graine"] for r in rows})
    methods = [m for m in METHODS if any(r["methode"] == m for r in rows)]
    idx = {(r["methode"], r["graine"]): r for r in rows}

    rep = {"repertoire": args.dir, "n_cellules": len(rows),
           "graines": seeds, "methodes": methods, "par_methode": {}}

    print(f"{args.dir}: {len(rows)} cellules, {len(seeds)} graines, {len(methods)} methodes\n")
    hdr = (f"{'methode':10s} {'steriles':>9s} {'faisable robuste':>22s} "
           f"{'marge et2 mediane [q1,q3]':>30s}")
    print(hdr); print("-" * len(hdr))

    marg_table, feas_table = {}, {}
    for m in methods:
        cells = [idx[(m, s)] for s in seeds if (m, s) in idx]
        n = len(cells)
        n_ster = sum(c["etage1"]["sterile"] for c in cells)
        n_feas = sum(bool(c["etage2"]["faisable"]) for c in cells)
        lo, hi = wilson(n_feas, n)
        marg = np.array([c["etage2"]["meilleure_marge"] for c in cells], float)
        marg_table[m] = np.array([idx[(m, s)]["etage2"]["meilleure_marge"]
                                  if (m, s) in idx else np.nan for s in seeds])
        feas_table[m] = np.array([float(idx[(m, s)]["etage2"]["faisable"])
                                  if (m, s) in idx else np.nan for s in seeds])
        q1, med, q3 = np.percentile(marg, [25, 50, 75])
        rep["par_methode"][m] = {
            "n": n, "steriles_etage1": n_ster,
            "faisables_robustes_conception": n_feas,
            "wilson95_faisable": [lo, hi],
            "marge_mediane": float(med), "marge_q1": float(q1), "marge_q3": float(q3),
            "actives": [c["etage2"].get("contrainte_active") for c in cells],
        }
        print(f"{m:10s} {n_ster:4d}/{n:<4d} {n_feas:3d}/{n:<3d}[{lo*100:5.1f},{hi*100:5.1f}]% "
              f"{med:+10.5f} [{q1:+.4f},{q3:+.4f}]")

    rep["statistiques_marge"] = compare_methods(marg_table, higher_is_better=True)
    st = rep["statistiques_marge"]
    print(f"\nGraines completes (unite statistique) : {st['n_graines_completes']}")
    if "friedman" in st and "p" in st["friedman"]:
        print(f"Friedman : chi2={st['friedman']['chi2']:.3f}  p={st['friedman']['p']:.3e}"
              f"  (Iman-Davenport F={st['friedman']['iman_davenport_F']:.3f})")
    sig = {k: v for k, v in st["paires"].items() if v["significatif_holm"]}
    print(f"Paires significatives apres Holm : {len(sig)} / {len(st['paires'])}")
    for k, v in sig.items():
        print(f"  {k:24s} p_holm={v['p_holm']:.3e}  effet={v['effet_rang_biserial']:+.3f}")

    out = args.out or os.path.join(args.dir, "_analyse.json")
    with open(out, "w") as f:
        json.dump(rep, f, indent=2)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
