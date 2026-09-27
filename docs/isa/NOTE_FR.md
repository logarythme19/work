# Note de synthèse (FR) — refonte ISA Transactions de CMO_LQI_PIDF_Buck_48V

Manuscrit : `docs/isa/main.tex` → `docs/isa/main.pdf` (elsarticle, anglais).
Toutes les valeurs du texte sont générées par `scripts/cmo_numbers.py` à partir de
`results/cmo/formula_first.json` et `results/cmo/analysis.json`.

## 1. Chaîne « formula first », de la Fig. 1 à la Fig. 2 et à l'équivalence exacte
| Étape | Résultat vérifié |
|---|---|
| Fig. 1, éq. (1)–(15) | équilibres, G_vd(s), écran CCM : identiques à la v6 (d_ff(1 V) = 0,029576, marge 2,96 mA, R_crit/R = 1,030) |
| Paramètres du `.slx` | décodés du *model workspace* (snubber Rs = 100 kΩ sur le MOSFET, diode de corps 0,01 Ω, PWM continue 50 kHz, Ts = 10 µs) |
| Fig. 2 (modèle commuté) | réplique Python : 15,0 % de DCM et 60,6 % de rendement à 1 V (Simulink v6 : 15,27 % et 60,51 %) ; bilan d'énergie < 5·10⁻⁶ % |
| LQI de Bryson | K = [0,0356434 ; 0,0121258 ; −98,34551], pôles −4243 et −7918 ± j6987 : identiques à la v6 |
| Proposition 1 (LQI → PIDF 2-DOF) | résidu de coefficients 0, erreur fréquentielle 8,7·10⁻¹⁶ |
| **Proposition 2 (nouvelle)** | l'état du filtre dérivé est exactement −v_C sur toute trajectoire : PIDF ≡ LQI sous saturation, anti-emballement, DCM, et pour tout L, r_L, R_on, r_d, V_f (écart 10⁻¹³ V). Écart explicite −K_ix(i_o − v_o/R) sinon |
| Échantillonnage | filtre discrétisé en BOZ : écart PIDF/LQI de 14,5 V sur l'échelon 1→46 V ; avec maintien triangulaire (FOH) : 0,35 V |
| Certificat LMI | α = 4157,93 rad/s, ζ = 0,73309 au point d'ancrage ; 6414,8 rad/s et 0,7496 pour le correcteur retenu |

## 2. Méthodologie à 7 algorithmes (PSO, GA, ABC, pAEABC, pIGWO, pIGWO-DLH + refineBO)
- 30 graines appariées, 240 + 60 appels, 33 inégalités en unités physiques dont 3 mesurées sur le modèle commuté.
- Les 6 branches convergent vers **le même optimum** (J_s = 0,8151 ; −18,5 % par rapport à Bryson).
- Elles diffèrent en fiabilité : pIGWO 28/30, pAEABC 25/30, GA 13/30 (Friedman p = 0,035).

## 3. Est-ce supérieur à PI, PID, PIDF, LQR, LQG, LQI sur 24 → 1 → 46 V ? (réplique commutée)
| Correcteur | Dépassement | i_L crête | Erreur max | Spécifications tenues |
|---|---|---|---|---|
| **CMO-LQI-PIDF (retenu, v_o seul)** | 0 V | 22,4 A | 25,6 mV (1 V) | 4/5 (manque de 0,6 mV) |
| CMO-LQI-SF (jumeau, i_L mesuré) | 0 V | 22,9 A | 0,0 mV | 5/5 |
| PI / PID / PIDF (règles standard) | 26,6–28,6 V | 40,8 A | ≤ 12 mV | 3/5 |
| LQR / LQG | 1,0 / 0,3 V | 24,4 A | 754 / 290 mV | 2/5, 3/5 |
| LQI (Bryson, état complet) | 9,4 V | 32,6 A | 0 mV | 2/5 |

**Réponse honnête :**
- **Oui** face aux six méthodes classiques : c'est la seule procédure qui supprime le dépassement,
  respecte le calibre de l'IRF540N et tient les niveaux ; les six comparateurs dépassent la
  spécification d'un facteur 6 à 30. Le suivi (IAE) reste à moins de 5 % du LQI.
- **Non** si l'on prétend que la paramétrisation LQI « optimise mieux » : la recherche directe
  des gains du même PIDF, avec le même protocole, atteint le même optimum (p = 0,74). L'apport
  de la paramétrisation LQI est structurel : chaque point est LQ-optimal (marges, certificat),
  alors que 42 % des solutions faisables de la recherche directe ne le sont pas.

## 4. À faire avant soumission
1. Rejouer `matlab/run_cmo_validation.m` sur le `.slx` et remplacer le Tableau 9 par les valeurs Simulink.
2. Compléter auteurs, affiliations, CRediT, et déposer le code (DOI).
3. Une validation expérimentale renforcerait fortement l'article pour ISA Transactions.
