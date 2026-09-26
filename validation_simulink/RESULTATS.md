# Résultats de la validation croisée Python ↔ Simulink

Exécutée par l'utilisateur sous MATLAB R2022a sur le modèle de référence
`buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx`, plant nominal, 6 candidats.
Tableau brut : `simulink_results.csv` ; comparaison : `comparaison.txt`.

## Verdict, selon les critères fixés avant l'essai
- **5 candidats sur 5 dotés d'un correcteur : ACCORD.** Écart sur J de 0,03 à 0,97 %,
  écart maximal sur les gᵢ ≤ 0,004, verdict de faisabilité identique. Les deux points
  faisables de la sonde restent faisables sous Simulink (+0,1911 et +0,0189).
- **Témoin sans rétroaction : DÉSACCORD sur g4** (|Δg4| = 0,0999 > 0,05). Il n'y a pas
  de changement de verdict (infaisable des deux côtés).

## Analyse du désaccord
L'écart porte sur l'erreur au palier de 1 V : e1 = 0,191 V (Python) contre 0,188 V
(Simulink), soit 3 mV, que la normalisation de g4 par 0,03 V amplifie en 0,10. Sans
rétroaction, rien ne corrige l'équilibre du convertisseur en conduction discontinue,
qui dépend du modèle de diode et d'interrupteur (loi affine d'un côté, blocs Simscape
avec diode de structure et snubber de l'autre). En boucle fermée, l'intégrateur efface
cet écart.

Même origine probable : Vmin est systématiquement plus haut de 2,3 mV sous Simulink
(0,7574 V contre 0,7551 V pour les meilleurs candidats des trois bras). **Le plancher de
tension à d_min reste sous 0,8 V dans Simulink** : le verrou physique du diagnostic E11
est confirmé par le modèle de référence.

## Portée
Plant nominal seulement ; coins d'incertitude non rejoués. Snubber supposé
(Rs = 1e5 Ω, Cs = ∞).
