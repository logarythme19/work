#!/usr/bin/env python3
"""Etage 3 (replay fin), validation, puis test SCELLE (par.5).

Ordre impose et irreversible:
  1. replay fin a 0.2 us des candidats retenus sur le jeu de CONCEPTION;
     aucun candidat n'est declare faisable sans ce replay effectif;
  2. evaluation sur le jeu de VALIDATION, qui n'a jamais servi a decider;
  3. test SCELLE: tirages uniformes, graine publiee avant evaluation, joue UNE
     SEULE FOIS en fin d'etude. Resultat = proportion + intervalle de Wilson,
     jamais "robuste".

Le test est protege par un fichier sentinelle: une seconde execution est
refusee, car rejouer le test apres avoir vu son resultat le transformerait en
jeu de conception (par.5).
"""

import argparse, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.evaluate import Controller, evaluate
from isafo.folqi import synthesize
from isafo.params import NUM
from isafo.problem import (corner_plants, random_plants, split_design_validation,
                           SEED_SPLIT, SEED_TEST)
from isafo.stats import wilson

LAB = ["g1", "g2", "g3", "g4", "g5", "g6"]


def controller_from_cell(cell):
    e2 = cell["etage2"]
    x = np.asarray(e2["x"], float)
    syn = synthesize(x[:7])
    theta = np.asarray(e2["theta"], float)
    return Controller.from_theta(theta, Kb=syn.Kb), theta


def run_set(ctrl, plants, h):
    res = [evaluate(ctrl, h=h, plant=P) for P in plants]
    k = sum(r.feasible for r in res)
    margins = [r.margin for r in res]
    active = [LAB[r.active_constraint() - 1] for r in res]
    return k, len(res), margins, active


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--n-test", type=int, default=400)
    ap.add_argument("--run-sealed-test", action="store_true",
                    help="joue le test scelle (UNE SEULE FOIS)")
    args = ap.parse_args()

    design, validation = split_design_validation(SEED_SPLIT)
    allp = corner_plants()
    cells = []
    for p in sorted(glob.glob(os.path.join(args.dir, "*.json"))):
        if os.path.basename(p).startswith("_"):
            continue
        with open(p) as f:
            c = json.load(f)
        if c["etage2"]["faisable"] and c["etage2"]["theta"] is not None:
            cells.append(c)

    print(f"{len(cells)} candidats faisables en conception a 1 us")
    out = {"candidats": []}
    for c in cells:
        ctrl, theta = controller_from_cell(c)
        kd_, nd_, md, ad = run_set(ctrl, [allp[i] for i in design], NUM.h_confirm)
        kv, nv, mv, av = run_set(ctrl, [allp[i] for i in validation], NUM.h_confirm)
        rec = {"methode": c["methode"], "graine": c["graine"], "bras": c.get("bras"),
               "theta": theta.tolist(),
               "conception_replay_fin": {"k": kd_, "n": nd_, "marge_min": min(md),
                                         "actives": ad},
               "validation_replay_fin": {"k": kv, "n": nv, "marge_min": min(mv),
                                         "actives": av,
                                         "wilson95": list(wilson(kv, nv))}}
        out["candidats"].append(rec)
        ok = "OK " if kd_ == nd_ else "KO "
        print(f"  {ok}{c['methode']:9s} {c['graine']}  conception {kd_}/{nd_} "
              f"(marge min {min(md):+.5f})   validation {kv}/{nv} "
              f"(marge min {min(mv):+.5f})")

    sentinel = os.path.join(args.dir, "_TEST_SCELLE_DEJA_JOUE")
    if args.run_sealed_test:
        if os.path.exists(sentinel):
            sys.exit("REFUS: le test scelle a deja ete joue. Le rejouer apres en "
                     "avoir vu le resultat en ferait un jeu de conception (par.5).")
        # Seul le candidat retenu AVANT le test y est soumis: le meilleur en
        # validation, departage par la marge de conception.
        ok_cands = [r for r in out["candidats"]
                    if r["conception_replay_fin"]["k"] == r["conception_replay_fin"]["n"]]
        if not ok_cands:
            print("\nAucun candidat faisable en conception apres replay fin:"
                  " le test scelle n'est PAS joue (il n'y a rien a tester).")
        else:
            best = max(ok_cands, key=lambda r: (r["validation_replay_fin"]["k"],
                                                r["conception_replay_fin"]["marge_min"]))
            cell = next(c for c in cells if c["methode"] == best["methode"]
                        and c["graine"] == best["graine"])
            ctrl, _ = controller_from_cell(cell)
            plants = random_plants(args.n_test, SEED_TEST)
            k, n, m, a = run_set(ctrl, plants, NUM.h_confirm)
            lo, hi = wilson(k, n)
            out["test_scelle"] = {"candidat": {"methode": best["methode"],
                                               "graine": best["graine"]},
                                  "graine": SEED_TEST, "k": k, "n": n,
                                  "wilson95": [lo, hi],
                                  "marge_min": min(m)}
            open(sentinel, "w").write(f"{SEED_TEST}\n")
            print(f"\nTEST SCELLE (graine {SEED_TEST}): {k}/{n} succes, "
                  f"Wilson 95 % = [{lo*100:.2f} % ; {hi*100:.2f} %]")

    with open(os.path.join(args.dir, "_replay_validation.json"), "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
