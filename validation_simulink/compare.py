#!/usr/bin/env python3
"""Confronte les chiffres Python (candidates.csv) et Simulink (simulink_results.csv).

Criteres d'accord FIXES AVANT d'avoir vu les resultats Simulink, pour ne pas
les ajuster apres coup :
  - verdict de faisabilite identique (signe de la marge) pour chaque candidat ;
  - ecart relatif sur J <= 5 % ;
  - ecart absolu sur chaque g_i <= 0,05 ;
  - ecart absolu sur Vmin <= 0,02 V et relatif sur Imax <= 5 %.
Les ecarts attendus viennent du solveur (Simscape contre RK4 a pas fixe), des
modeles de semi-conducteurs et du snubber, absent du modele Python.
"""
import csv, os, sys

here = os.path.dirname(os.path.abspath(__file__))
py = {r["nom"]: r for r in csv.DictReader(open(os.path.join(here, "candidates.csv")))}
p_sl = os.path.join(here, "simulink_results.csv")
if not os.path.exists(p_sl):
    sys.exit("simulink_results.csv absent : lancer d'abord run_crossval.m sous MATLAB.")
sl = {r["nom"]: r for r in csv.DictReader(open(p_sl))}

f = lambda r, k: float(r[k])
ok_all = True
print(f"{'candidat':18s} {'J py':>8s} {'J sl':>8s} {'dJ %':>6s} {'m py':>8s} {'m sl':>8s} "
      f"{'max|dg|':>8s} {'Vmin py/sl':>15s} {'verdict':>9s}")
for nom, a in py.items():
    if nom not in sl:
        print(f"{nom:18s} absent des resultats Simulink"); ok_all = False; continue
    b = sl[nom]
    gp = [f(a, f"g{i}_py") for i in range(1, 7)]
    gs = [f(b, f"g{i}_sl") for i in range(1, 7)]
    mp, ms = -max(gp), -max(gs)
    dJ = abs(f(b, "J_sl") - f(a, "J_py")) / abs(f(a, "J_py"))
    dg = max(abs(x - y) for x, y in zip(gp, gs))
    dV = abs(f(b, "Vmin_sl") - f(a, "Vmin_py"))
    dI = abs(f(b, "Imax_sl") - f(a, "Imax_py")) / f(a, "Imax_py")
    same = (mp >= 0) == (ms >= 0)
    ok = same and dJ <= 0.05 and dg <= 0.05 and dV <= 0.02 and dI <= 0.05
    ok_all &= ok
    print(f"{nom:18s} {f(a,'J_py'):8.5f} {f(b,'J_sl'):8.5f} {dJ*100:6.2f} {mp:+8.4f} {ms:+8.4f} "
          f"{dg:8.4f} {f(a,'Vmin_py'):7.4f}/{f(b,'Vmin_sl'):7.4f} "
          f"{('ACCORD' if ok else ('DESACCORD' if same else 'VERDICT!')):>9s}")
print("\nConclusion :", "accord sur tous les candidats" if ok_all else
      "au moins un desaccord : a analyser avant toute conclusion")
