#!/usr/bin/env python3
"""Exporte les candidats a rejouer sous Simulink, avec leurs chiffres Python.

Chaque candidat est evalue ici au pas de confirmation (0,2 us) sur le plant
NOMINAL, avec exactement le meme calcul de J et des six contraintes que la
campagne. Le script MATLAB rejoue les memes candidats sur le modele Simulink ;
validation_simulink/compare.py confronte ensuite les deux.
"""
import csv, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from isafo.evaluate import Controller, evaluate
from isafo.folqi import synthesize
from isafo.params import NUM

rows = []
p = json.load(open("results/probe_setpoint_weight.json"))
for key, name in (("b_libre", "sonde_b_libre"), ("b_fixe_1_proposition_1", "sonde_b1")):
    th = np.array(p["bras"][key]["theta_meilleur"])
    rows.append((name, th[:7], 1.0 / 10 ** th[7]))
for arm in ("v5", "v5b", "v5kp"):
    cells = [json.load(open(f)) for f in glob.glob(f"results/campaign_{arm}/*.json")
             if not os.path.basename(f).startswith("_")]
    best = max(cells, key=lambda c: c["etage2"]["meilleure_marge"])
    rows.append((f"meilleur_{arm}", np.array(best["etage2"]["theta"]),
                 synthesize(np.array(best["etage2"]["x"][:7])).Kb))
rows.append(("feedforward_seul", np.array([-9, -9, -9, 1, 1, np.log10(2e4), 1]), 1e3))

out = "validation_simulink/candidates.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["nom", "kp", "ki", "kd", "lambda", "mu", "wf", "b", "c", "Kb",
                "J_py", "g1_py", "g2_py", "g3_py", "g4_py", "g5_py", "g6_py",
                "Imax_py", "Vmax_py", "Vmin_py", "e24_py", "e1_py", "e46_py"])
    for name, th, Kb in rows:
        c = Controller.from_theta(th, Kb=Kb)
        if name == "feedforward_seul":
            c.kd = 0.0
        r = evaluate(c, h=NUM.h_confirm)
        w.writerow([name, c.kp, c.ki, c.kd, c.lam, c.mu, c.wf, c.b, c.c, c.Kb, r.J,
                    *r.g.tolist(), r.Imax, r.Vmax, r.Vmin, r.e24, r.e1, r.e46])
        print(f"{name:18s} J={r.J:.5f} marge={r.margin:+.4f} Imax={r.Imax:.3f} Vmin={r.Vmin:.4f}")
print(f"-> {out}")
