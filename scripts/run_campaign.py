#!/usr/bin/env python3
"""Recherche etagee, six methodes, graines appariees (par.3.3).

Etage 1 : faisabilite nominale, un seul scenario, 700 appels. Critere de
          passage: au moins 20 points faisables DISTINCTS, sinon la graine est
          declaree STERILE et rapportee comme telle (jamais silencieusement
          ecartee: le par.6 impose de conserver tous les echecs).
Etage 2 : robustesse sur le jeu de CONCEPTION, 300 appels, amorces sur les
          points faisables de l'etage 1 DE LA MEME GRAINE, jamais sur un
          candidat archive d'une campagne anterieure.
Etage 3 : replay fin a 0.2 us (script separe), car aucun candidat n'est declare
          faisable sans replay fin effectif (par.5).
"""

import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.algos import run, METHODS, DEFAULT_POP
from isafo.anchor import SaturationAudit
from isafo.deb import best_of
from isafo.folqi import V5_LB, V5_UB
from isafo.params import NUM
from isafo.problem import (Problem, corner_plants, split_design_validation,
                           SEED_BASE, SEED_SPLIT, V5B_LB, V5B_UB,
                           V5KP_LB, V5KP_UB)

BUDGET1, BUDGET2, MIN_FEASIBLE = 700, 300, 20


def _distinct(xs, tol=4):
    return len({tuple(np.round(x, tol)) for x in xs})


def one_cell(method: str, seed: int, budget1: int, budget2: int,
             design: list, arm: str = "v5") -> dict:
    """Une cellule (methode, graine): etage 1 puis etage 2."""
    t0 = time.time()
    free_b = (arm == "v5b")
    free_kp = (arm == "v5kp")
    LB, UB = {"v5": (V5_LB, V5_UB), "v5b": (V5B_LB, V5B_UB),
              "v5kp": (V5KP_LB, V5KP_UB)}[arm]
    audit1 = SaturationAudit()
    prob1 = Problem(cases=[corner_plants()[0]], audit=audit1, free_b=free_b,
                    free_kp=free_kp)

    h1 = run(method, prob1, budget1, seed, LB, UB,
             n_pop=DEFAULT_POP, mode="margin", stage="etage1")

    feas1 = [r for r in h1 if r.feasible]
    n_distinct = _distinct([r.x for r in feas1])
    sterile = n_distinct < MIN_FEASIBLE

    b1 = best_of(h1, "margin")
    out = {
        "methode": method, "graine": seed, "bras": arm,
        "etage1": {
            "appels": len(h1),
            "n_faisables": len(feas1),
            "n_faisables_distincts": n_distinct,
            "sterile": bool(sterile),
            "meilleure_marge": float(b1.margin) if b1 else None,
            "meilleur_J": float(b1.J) if b1 else None,
            "meilleure_violation": float(b1.violation) if b1 else None,
            "x": b1.x.tolist() if b1 else None,
            "saturation": audit1.to_dict(),
        },
    }

    # Amorcage de l'etage 2 sur les faisables de l'etage 1 de la MEME graine.
    order = np.argsort([-r.margin for r in feas1])
    seeds2 = np.array([feas1[i].x for i in order[:DEFAULT_POP]]) if feas1 else None

    audit2 = SaturationAudit()
    prob2 = Problem(cases=[corner_plants()[i] for i in design], audit=audit2,
                    free_b=free_b, free_kp=free_kp)
    h2 = run(method, prob2, budget2, seed + 1, LB, UB,
             n_pop=DEFAULT_POP, seeds=seeds2, mode="margin", stage="etage2")

    feas2 = [r for r in h2 if r.feasible]
    b2 = best_of(h2, "margin")
    out["etage2"] = {
        "appels": len(h2),
        "n_faisables": len(feas2),
        "amorce_sur": int(len(seeds2)) if seeds2 is not None else 0,
        "meilleure_marge": float(b2.margin) if b2 else None,
        "meilleur_J": float(b2.J) if b2 else None,
        "meilleure_violation": float(b2.violation) if b2 else None,
        "x": b2.x.tolist() if b2 else None,
        "theta": b2.theta.tolist() if (b2 and b2.theta is not None) else None,
        "faisable": bool(b2.feasible) if b2 else False,
        "contrainte_active": int(b2.active_constraint()) if b2 else None,
        "saturation": audit2.to_dict(),
    }
    if b2 is not None and "par_cas" in b2.extra:
        out["etage2"]["actives_par_cas"] = [c["active"] for c in b2.extra["par_cas"]]
        out["etage2"]["marges_par_cas"] = [c["margin"] for c in b2.extra["par_cas"]]
    out["secondes"] = time.time() - t0
    return out


def _job(a):
    return one_cell(*a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=30)
    ap.add_argument("--methods", nargs="*", default=list(METHODS))
    ap.add_argument("--budget1", type=int, default=BUDGET1)
    ap.add_argument("--budget2", type=int, default=BUDGET2)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--arm", choices=("v5", "v5b", "v5kp"), default="v5")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    design, _ = split_design_validation(SEED_SPLIT)
    args.out = args.out or f"results/campaign_{args.arm}"
    os.makedirs(args.out, exist_ok=True)

    jobs = []
    for r in range(args.seeds):
        for m in args.methods:
            seed = SEED_BASE + r
            p = os.path.join(args.out, f"{m}_{seed}.json")
            if not os.path.exists(p):
                jobs.append((m, seed, args.budget1, args.budget2, design, args.arm))

    print(f"bras {args.arm}: {len(jobs)} cellules a calculer sur {args.workers} coeurs")
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_job, j): j for j in jobs}
        for fut in as_completed(futs):
            m, seed = futs[fut][0], futs[fut][1]
            try:
                res = fut.result()
            except Exception as exc:
                print(f"  ECHEC {m} {seed}: {exc}")
                continue
            with open(os.path.join(args.out, f"{m}_{seed}.json"), "w") as f:
                json.dump(res, f, indent=2)
            done += 1
            e1, e2 = res["etage1"], res["etage2"]
            tag = "STERILE" if e1["sterile"] else f"{e1['n_faisables_distincts']:3d} faisables"
            el = time.time() - t0
            print(f"[{done:4d}/{len(jobs)}] {m:9s} {seed}  et1: {tag:14s}"
                  f"  et2: marge={e2['meilleure_marge']:+.5f}"
                  f" J={e2['meilleur_J'] if e2['meilleur_J'] else float('nan'):.4f}"
                  f"  ({res['secondes']:.0f}s, ecoule {el/60:.1f}min)")

    print(f"\nTermine en {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
