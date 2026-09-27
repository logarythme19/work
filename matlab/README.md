# Rejouer les résultats sur le modèle Simulink de référence

Fichiers :
- `cmo_sfun.m` — S-function niveau 2, portage ligne à ligne de `cmo/ctrl.py`
  (PI, PID, PIDF 2-DOF avec filtre à maintien triangulaire, LQR, LQG, LQI).
  Même interface que `buck_comparator_sfun` : entrées Vref, V_o, I_L, I_c, V_in ;
  sortie d_sat avec retard d'un échantillon.
- `run_cmo_validation.m` — rejoue chaque ligne de `results/cmo/controllers.csv`
  sur le `.slx` sans jamais le modifier (Simulink.SimulationInput), calcule IAE,
  dépassement, courant crête, erreurs statiques, rendement, bilan d'énergie.
- `cmo_formula_first.m` — vérification indépendante de la chaîne analytique
  (équilibre, G_vd, CARE, gain LQI, correspondance exacte LQI → PIDF 2-DOF).

Utilisation (R2022a+, Simscape Electrical, Control System Toolbox) :
```matlab
cd <repo>/matlab
cmo_formula_first
T = run_cmo_validation('buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx', ...
                       '../results/cmo/controllers.csv');
```
Le résultat est écrit dans `matlab/cmo_validation_simulink.csv`. Comparer avec
`results/cmo/analysis.json` (clé `switched`). Conditions initiales imposées :
iL(0) = 2,4 A, vC(0) = 24 V (le `.slx` sauvegardé part de 1,2 A / 12 V).
