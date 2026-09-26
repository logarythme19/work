#!/usr/bin/env python3
"""Mesure de l'impact du defaut de feedforward (journal E12).

La campagne v6.1 a calcule le feedforward avec les parametres du plant simule
au lieu des valeurs nominales figees (par.1). Seul l'etage 2 (coins
d'incertitude) est concerne ; l'etage 1 ne simule que le plant nominal.

On reevalue ici, avec le noyau corrige, le candidat retenu de chaque cellule
sur le jeu de conception, et on compare a la marge enregistree. Limite : la
RECHERCHE de l'etage 2 aurait pu suivre une autre trajectoire avec le noyau
corrige ; cette mesure ne porte que sur les candidats retenus.
"""
import glob, json, os, sys
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np


def one(f):
    from isafo.problem import Problem, corner_plants, split_design_validation
    c = json.load(open(f))
    arm = c.get("bras", "v5")
    x = np.array(c["etage2"]["x"])
    prob = Problem(cases=[corner_plants()[i] for i in split_design_validation()[0]],
                   free_b=(arm == "v5b"), free_kp=(arm == "v5kp"))
    r = prob(x)
    return arm, c["etage2"]["meilleure_marge"], r.margin, r.feasible


files = [f for a in ("v5", "v5b", "v5kp") for f in glob.glob(f"results/campaign_{a}/*.json")
         if not os.path.basename(f).startswith("_")]
with ProcessPoolExecutor(4) as ex:
    res = list(ex.map(one, files, chunksize=8))
out = {}
for arm in ("v5", "v5b", "v5kp"):
    R = [(a, b, c, d) for a, b, c, d in res if a == arm]
    d = np.array([c - b for _, b, c, _ in R])
    out[arm] = {"n": len(R), "delta_marge_max_abs": float(np.max(np.abs(d))),
                "delta_marge_median": float(np.median(d)),
                "faisables_apres_correction": int(sum(x[3] for x in R)),
                "meilleure_marge_avant": float(max(x[1] for x in R)),
                "meilleure_marge_apres": float(max(x[2] for x in R))}
    o = out[arm]
    print(f"{arm:5s} n={o['n']}  |delta marge| max={o['delta_marge_max_abs']:.2e}  mediane={o['delta_marge_median']:+.2e}"
          f"  faisables apres correction={o['faisables_apres_correction']}"
          f"  meilleure marge {o['meilleure_marge_avant']:+.4f} -> {o['meilleure_marge_apres']:+.4f}")
json.dump(out, open("results/impact_feedforward.json", "w"), indent=2)
