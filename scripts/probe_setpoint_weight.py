#!/usr/bin/env python3
"""Sonde decisive: b est-il le levier qui separe faisable et infaisable ?

La proposition 1 fixe b = 1 EXACTEMENT, et le vecteur v5 n'a aucune composante
agissant sur b. Si le meilleur FO-PIDF avec b = 1 reste infaisable alors qu'il
existe des FO-PIDF faisables avec b libre, la parametrisation v5 exclut
structurellement la region faisable, independamment de la qualite de la
recherche.

Ce n'est PAS un balayage mono-varie (par.4.3): chaque bras est une recherche
CONJOINTE sur tout theta, seule la valeur de b differe entre les bras.
"""

import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.algos import run
from isafo.deb import Record
from isafo.evaluate import Controller, evaluate
from isafo.params import NUM

LB = np.array([-2.0, -1.0, -6.0, 0.40, 0.40, np.log10(500), 0.0, -4.3])
UB = np.array([2.5, 5.0, -1.0, 1.00, 1.00, np.log10(2e4), 1.0, -2.3])
LAB = ["g1 Imax", "g2 Vmax", "g3 e24", "g4 e1", "g5 e46", "g6 Vmin"]


def make_fn(b_fixed):
    def fn(x):
        x = np.array(x, float)
        if b_fixed is not None:
            x[6] = b_fixed
        c = Controller.from_theta(x[:7], Kb=1.0 / 10 ** x[7])
        r = evaluate(c)
        return Record(x=x, J=r.J, g=r.g, margin=r.margin, feasible=r.feasible,
                      violation=r.violation,
                      extra={"Imax": r.Imax, "Vmax": r.Vmax, "Vmin": r.Vmin})
    return fn


def arm(label, b_fixed, budget, methods_seeds):
    best = None
    rows = []
    for m, s in methods_seeds:
        h = run(m, make_fn(b_fixed), budget, s, LB, UB, n_pop=12, mode="margin")
        b = max(h, key=lambda r: (r.feasible, r.margin))
        n_feas = sum(r.feasible for r in h)
        rows.append({"methode": m, "graine": s, "marge": b.margin,
                     "n_faisables": n_feas, "Imax": b.extra["Imax"],
                     "b": float(b.x[6]), "x": b.x.tolist(),
                     "active": LAB[int(np.argmax(b.g))]})
        if best is None or (b.feasible, b.margin) > (best.feasible, best.margin):
            best = b
    return best, rows


def main():
    ms = [("IGWO_DLH", 1), ("PSO", 2), ("AEABC", 3), ("GA", 4),
          ("IGWO", 5), ("ABC", 6)]
    budget = 1500
    out = {"budget_par_lancer": budget, "bras": {}}

    for label, bf in (("b_libre", None), ("b_fixe_1_proposition_1", 1.0)):
        t0 = time.time()
        best, rows = arm(label, bf, budget, ms)
        print(f"\n=== {label}  ({time.time()-t0:.0f}s) ===")
        for r in rows:
            print(f"  {r['methode']:9s} marge={r['marge']:+.5f}  faisables={r['n_faisables']:4d}"
                  f"  Imax={r['Imax']:6.2f}  b={r['b']:.3f}  active={r['active']}")
        print(f"  -> MEILLEURE MARGE = {best.margin:+.5f}")

        # Replay fin (par.5): aucun candidat faisable sans replay a 0.2 us.
        c = Controller.from_theta(best.x[:7], Kb=1.0 / 10 ** best.x[7])
        rf = evaluate(c, h=NUM.h_confirm)
        print(f"  replay fin 0.2 us : marge={rf.margin:+.5f}  faisable={rf.feasible}"
              f"  (ecart de marge {abs(rf.margin-best.margin):.2e})")

        out["bras"][label] = {
            "lancers": rows,
            "meilleure_marge_1us": best.margin,
            "meilleure_marge_02us": rf.margin,
            "faisable_apres_replay_fin": bool(rf.feasible),
            "theta_meilleur": best.x.tolist(),
            "g_replay_fin": rf.g.tolist(),
        }

    os.makedirs("results", exist_ok=True)
    with open("results/probe_setpoint_weight.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\n-> results/probe_setpoint_weight.json")


if __name__ == "__main__":
    main()
