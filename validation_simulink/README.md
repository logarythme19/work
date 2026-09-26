# Validation croisée Python ↔ Simulink

Rejoue sur le modèle Simulink de référence les candidats figés de la campagne,
avec **le même correcteur** et **les mêmes formules** de coût et de contraintes,
puis compare aux chiffres Python.

## Fichiers
- `candidates.csv` : 6 candidats (2 faisables issus de la sonde, le meilleur de chaque
  bras, un témoin sans rétroaction), avec leurs chiffres Python au pas de 0,2 µs.
  Il est régénéré par `python scripts/export_candidates.py`.
- `fopidf_sfun.m` : S-function niveau 2, traduction ligne à ligne de `isafo/_kernel.py`,
  avec la même interface que `buck_comparator_sfun`.
- `run_crossval.m` : rejoue les candidats sans jamais modifier le `.slx`.
- `compare.py` : confronte les deux jeux de chiffres selon des critères fixés à l'avance.

## Utilisation
Sous MATLAB (R2022a ou plus récent, avec Simscape Electrical) :
```matlab
cd ~/work/validation_simulink
T = run_crossval('/chemin/vers/buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx');
```
Puis, dans le terminal :
```bash
cd ~/work && python validation_simulink/compare.py
```

## Hypothèses à vérifier
- **Snubber** : Rs = 1e5 Ω et Cs = inf (snubber résistif), car les valeurs ne figurent pas
  dans le `.slx`. Remplacez-les dans `run_crossval.m` si vous les connaissez.
- **Plant nominal seulement.** Les coins d'incertitude ne sont pas rejoués.
- Le terme en |Δd| de J est calculé sur le rapport cyclique **appliqué** (retardé d'un
  échantillon) ; son poids (0,001) rend l'écart négligeable.
