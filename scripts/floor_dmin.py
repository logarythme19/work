#!/usr/bin/env python3
"""Diagnostic POST HOC (journal E11) : plancher de tension a d = d_min.

Constat de campagne : les marges g6 reviennent a l'identique d'un correcteur a
l'autre (Vmin = 0,7541 V au nominal). Hypothese testee ici : c'est l'equilibre
en conduction discontinue du convertisseur a rapport cyclique fixe d_min. Si
la commande reste saturee au plancher pendant la descente 24 -> 1 V, la sortie
s'y stabilise, quel que soit le correcteur.

Rapporte, jamais utilise pour selectionner (par.8).
"""
import json, os, sys
from dataclasses import replace
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from isafo.params import MISSION, NUM, PLANT
from isafo.evaluate import Controller, evaluate
from isafo.problem import corner_plants, split_design_validation

out = {"d_min": NUM.duty_min, "seuil_g6_V": 0.8, "cas": {}}
design, validation = split_design_validation()
for i, Pk in enumerate(corner_plants()):
    # Le feedforward est nominal fige (par.1, E12) : la consigne qui donne
    # exactement d = d_min se calcule avec le plant NOMINAL.
    P0 = PLANT
    r = NUM.duty_min * P0.R * P0.Vin / ((P0.R + P0.rL + P0.rd) + NUM.duty_min * (P0.Ron - P0.rd))
    M = replace(MISSION, v0=r, v1=r, v2=r, ip_amp=0.0, iL0=0.0, vC0=0.9)
    res = evaluate(Controller(kp=1e-12, ki=1e-12, kd=0.0, lam=1, mu=1, wf=2e4, b=1, Kb=0.0),
                   h=NUM.h_confirm, plant=Pk, mission=M)
    jeu = "conception" if i in design else "validation"
    out["cas"][i] = {"jeu": jeu, "R": Pk.R, "L": Pk.L, "C": Pk.C, "rC": Pk.rC,
                     "Vmin_plancher_V": res.Vmin, "g6_au_plancher": (0.8 - res.Vmin) / 0.8}
    print(f"cas {i:2d} ({jeu:10s}) R={Pk.R:5.1f}  plancher Vmin={res.Vmin:.4f} V  g6={(0.8-res.Vmin)/0.8:+.4f}")
n = sum(v["g6_au_plancher"] > 0 for v in out["cas"].values())
out["n_cas_plancher_sous_seuil"] = n
print(f"\n{n}/{len(out['cas'])} cas ont un plancher sous 0,8 V")
json.dump(out, open("results/plancher_dmin.json", "w"), indent=2)
