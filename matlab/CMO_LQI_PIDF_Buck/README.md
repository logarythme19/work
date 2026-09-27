# CMO-LQI-PIDF — Convertisseur Buck non idéal : fichiers MATLAB / Simulink

Ces fichiers accompagnent l'article « Structural LQI–PIDF equivalence and constrained
weight tuning for a wide-range Buck converter ». Ils permettent :

1. de **tester les paramètres optimaux** publiés sur le modèle Simulink de référence
   et sur une réplique commutée écrite en MATLAB ;
2. de **relancer le protocole de recherche** (6 métaheuristiques + refineBO) ;
3. de faire tout cela **pour un autre convertisseur** (Vin, L, rL, C, rC, R, Ron,
   rd, Vf, fs, Ts, niveaux de consigne, cahier des charges...).

## Contenu

| Fichier | Rôle |
|---|---|
| `cmo_params.m` | **Tous les paramètres** : étage de puissance, MLI, niveaux, mission, cahier des charges (33 inégalités), synthèse LQI, poids optimaux publiés, réglages de la recherche, noms des blocs Simulink. |
| `cmo_lib.m` | **Toutes les fonctions** : modèle non idéal, CARE/DARE, LQI → PIDF 2-DOF exact, comparateurs PI/PID/PIDF/LQR/LQG/LQI, couche de commande, modèle moyen, réplique commutée, évaluateur (6 scénarios, J_s, 33 contraintes), algorithmes 1–7, rejeu Simulink. |
| `cmo_sfun.m` | **S-function MATLAB** (niveau 2) : couche de commande échantillonnée commune, avec la même interface que le bloc « Common sampled controller » du modèle. |
| `buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx` | **Modèle Simulink de référence** (R2022a, Specialized Power Systems). Il n'est jamais modifié : les paramètres sont injectés au moment de la simulation. |
| `run_optimal.m` | Teste les paramètres optimaux et les 6 comparateurs (réplique + Simulink). |
| `run_search.m` | Lance le protocole de recherche des paramètres optimaux. |
| `run_formula_first.m` | Chaîne analytique : équilibres, CCM/DCM, Gvd(s), CARE, réalisation PIDF, identités. |

## Prérequis

- MATLAB R2019b ou plus récent. Les fonctions de calcul n'utilisent **aucune boîte à outils**.
- Pour le rejeu sur le `.slx` : Simulink et Simscape Electrical (Specialized Power Systems), R2022a ou plus récent.
- Facultatif : Parallel Computing Toolbox (`use_parfor = true` dans `run_search.m`).

Placez tous les fichiers dans **un même dossier**, puis ouvrez ce dossier dans MATLAB.

## 1. Tester les paramètres optimaux

```matlab
run_formula_first      % vérifications analytiques (quelques secondes)
run_optimal            % réplique commutée + modèle Simulink
```

`run_optimal` affiche le gain LQI et les gains du PIDF 2-DOF équivalent. Pour le convertisseur
nominal, on obtient Kp = 0,17436, Ki = 945,81, Kd = 8,2629e-6, N = 2e5 rad/s, b = 1, c = 0
et Kb = 92 951. Le script rejoue ensuite CMO-LQI-PIDF, son jumeau à retour d'état
(CMO-LQI-SF) et les comparateurs PI, PID, PIDF, LQR, LQG et LQI sur la mission
24 → 1 → 46 V avec impulsion de charge. Il écrit :

- `cmo_controllers.csv` : les vecteurs de paramètres `p` de chaque correcteur ;
- `cmo_validation_simulink.csv` : les indices mesurés sur le `.slx` (IAE, courant crête, dépassement, erreurs statiques, rendement, résidu du bilan d'énergie).

Pour observer le modèle dans Simulink, mettez `open_model = true` dans `run_optimal.m`.
Le modèle s'ouvre alors configuré en mémoire (le fichier n'est pas sauvegardé) : cliquez
sur **Run** et ouvrez les Scopes.

> N'exécutez pas le `.slx` directement sans ce script : il appelle une S-function
> d'origine (`buck_comparator_sfun`) qui n'est pas fournie. `run_optimal`,
> `F.sim_simulink` et `F.configure_model` la remplacent par `cmo_sfun`.

## 2. Relancer la recherche des paramètres optimaux

Ouvrez `run_search.m`, puis choisissez :

- `family` : `CMO-LQI-PIDF` (proposé), `CMO-LQI-SF`, `PIDF-direct`, `PIDF7-direct`, `PID-direct` ou `PI-direct` ;
- `methods` : une partie ou la totalité de `{'PSO','GA','ABC','pAEABC','pIGWO','pIGWO-DLH'}` ;
- `seeds` : les graines. L'article en utilise 30 : `P.search.seeds` ;
- le budget : `P.search.B_global = 240` et `P.search.B_local = 60` (valeurs de l'article).

Chaque appel simule 6 scénarios moyennés et un écran commuté de 34 ms, soit quelques
secondes sous MATLAB. Une exécution de 300 appels prend de l'ordre de 10 à 20 minutes.
Le résultat est enregistré dans `cmo_search_result.mat` (variable `best`). Pour rejouer
le correcteur trouvé :

```matlab
load cmo_search_result.mat best
% dans run_optimal.m :  xi = best.x;
```

Les générateurs aléatoires de MATLAB et de NumPy diffèrent. Les trajectoires de recherche
ne sont donc pas identiques bit à bit à celles de l'article, mais le protocole est le
même : populations, opérateurs, budget, règle de Deb, refineBO. L'évaluateur, lui,
reproduit exactement l'article : au point optimal publié, il donne J_s = 0,815095 et les
33 contraintes aux mêmes valeurs à 8 chiffres près.

## 3. Étudier un autre convertisseur

Il suffit de modifier `cmo_params.m`, ou de passer des paires nom/valeur à la place de la
ligne `P = cmo_params();` dans les scripts :

```matlab
P = cmo_params('Vin', 36, 'L', 47e-6, 'C', 220e-6, 'rC', 0.02, 'R', 5, ...
               'Vmid', 18, 'Vhigh', 34, 'spec.I_LIM', 15);
```

Les grandeurs suivantes sont recalculées automatiquement à partir des valeurs choisies :

- l'ancre de Bryson xi0 ;
- les bornes de recherche (log10(π/Ts)) ;
- la demi-ondulation de courant ;
- les six scénarios, l'écran commuté et la mission Simulink.

Sur le `.slx`, les blocs L, C, R, rl, R1, IGBT/Diode, Diode1, la source Vin et la
fréquence de la MLI reçoivent les nouvelles valeurs, de même que les conditions initiales
(i_L = Vmid/R, v_C = Vmid) et les profils de consigne et de charge.

Les poids optimaux publiés (`P.xi_opt`) valent pour le convertisseur nominal. Pour un autre
convertisseur, lancez `run_search` et utilisez le résultat obtenu.

Conditions de validité de l'équivalence exacte (Proposition 2 de l'article) :

- le réseau de sortie (R, rC, C) est connu ;
- N = 1/(rC·C) < π/Ts ;
- le PIDF ne lit que v_o. L'écart résiduel avec le LQI à état complet vaut −Ki_x·(i_o − v_o/R) et n'apparaît que sous perturbation de charge.

## Vecteur de paramètres du correcteur `p` (36 valeurs)

| Indice | Contenu |
|---|---|
| p(1) | code : 1 PI, 2 PID, 3 PIDF 2-DOF (CMO-LQI-PIDF), 4 LQR, 5 LQG, 6 LQI |
| p(2:8) | PI : Kp Ki Kb · PID : Kp Ki Kd b Kb · PIDF : Kp Ki Kd N b c Kb |
| p(9:20) | gains LQR/LQG/LQI, observateur de Kalman (LQG) |
| p(25:33) | R rL rd Ron Vin Vf rC C Ts nominaux (anticipation, jamais le procédé) |
| p(34:35) | dmin dmax |
| p(36) | niveau de consigne initial (initialisation de la S-function) |
