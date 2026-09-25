# Journal des écarts au cadrage

Le §0 du cadrage impose que tout écart soit **déclaré et justifié, jamais silencieux**.
Ce journal les recense tous, dans l'ordre chronologique. Les amendements de
protocole (A1–A4) sont aussi consignés dans `protocol/protocol_v6_1.json`.

| # | Écart | Nature | Justification | Moment |
|---|---|---|---|---|
| E1 | Chaîne reconstruite en Python, sans `anchor_projection.m` | outil | Code antérieur inaccessible dans la session ; choix validé par l'utilisateur. La projection est re-dérivée de la proposition 1 (article antérieur) et vérifiée à 5,4·10⁻¹⁶, au lieu de réutiliser un ancrage empirique figé indisponible. | avant tout calcul |
| E2 | Évaluateur **commuté** (diode bloquante, DCM) plutôt que moyenné | modèle numérique | L'espace de travail du `.slx` documente ~3,2 ms de conduction discontinue au saut 24 → 1 V. Le **plant** reste celui du §1, d'ordre entier. | avant tout calcul |
| E3 | Découpage exact des fronts MLI dans l'intégrateur | modèle numérique | Sinon 5 % de résolution en rapport cyclique à 1 µs. Conséquence : l'écart 1 µs → 0,2 µs n'est **pas** comparable aux −21 %/−12 %/+0,08 % de la campagne antérieure. | avant tout calcul |
| E4 | Intégrale réalisée comme s⁻¹·s^(1−λ) | réalisation | Seule façon d'obtenir une erreur statique nulle avec une bande Oustaloup commençant à 1 rad/s ; λ = 1 redonne l'intégrateur exact. | avant tout calcul |
| E5 | Ensemble d'incertitude choisi (L, C, R ±20 %, rC ±50 %) | protocole | Le cadrage cite 17 cas / 16 coins sans les définir. | protocole v6 |
| E6 | Bornes hautes d'ancrage non citées par le cadrage fixées non serrantes | protocole | Seules lb et ub(kp) = +2 sont données au §2. | protocole v6 |
| E7 | Jeu de test 300 → 400 | protocole (A1) | Wilson exact : 300/300 donne 98,74 %, 381 sont requis pour 99 %. | protocole v6.1 |
| E8 | Bras v5b (b libre) et v5kp (kp libre) | protocole (A2, A3) | v5 conforme reste le bras principal ; les deux autres sont des extensions/ablations déclarées. | protocole v6.1 |
| E9 | Snubbers du `.slx` non reproduits | modèle | Valeurs Rs, Cs absentes du fichier ; n'interviennent ni dans J ni dans g. | déclaré en limite |
| E10 | **Diagnostic d'oscillation ajouté après observation** | analyse | Cycle limite (7,7 kHz, 624 mV cc au palier 24 V) observé sur un candidat faisable. **Rapporté, jamais utilisé pour sélectionner** (§8). | après la sonde, pendant la campagne v6.1 |

## Hypothèses émises puis réfutées

- *« b = 1 (proposition 1) exclut la région faisable. »* — **Réfutée** par la sonde :
  avec b = 1, une recherche conjointe sur θ trouve des FO-PIDF faisables (marge +0,016
  après replay fin). b est un levier de marge, pas une barrière.
- *« Le cycle du rapport cyclique est une bascule à la fréquence de Nyquist. »* — **Partiellement
  réfutée** : le candidat b = 1 montre une alternance d'échantillon à échantillon de faible amplitude,
  mais le candidat b libre présente un cycle limite à 7,7 kHz, bien en dessous de Nyquist.

## Diagnostics post hoc ajoutés après la campagne v6.1

| # | Diagnostic | Constat | Statut |
|---|---|---|---|
| E11 | Plancher de tension à d = d_min (`scripts/floor_dmin.py`) | Les marges g6 reviennent à l'identique d'un correcteur à l'autre (Vmin = 0,7541 V au nominal). En boucle ouverte à d = 0,02, la sortie se stabilise sous 0,8 V pour la plupart des cas : le seuil de g6 est **au-dessus** du plancher physique du convertisseur saturé. | Rapporté, jamais utilisé pour sélectionner (§8). |
