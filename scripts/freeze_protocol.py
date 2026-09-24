#!/usr/bin/env python3
"""Livrable 1: protocole fige en JSON, avec empreintes SHA-256, AVANT calcul.

Le par.9 impose cet ordre et le par.8 interdit de modifier une empreinte pour
dissimuler un changement de code. Le protocole embarque donc l'empreinte de
chaque fichier source de la chaine d'evaluation: toute modification ulterieure
d'un de ces fichiers rend l'empreinte du protocole caduque et VISIBLE.
"""

import hashlib, json, os, sys, datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.params import invariants_dict, NUM
from isafo.folqi import V5_NAMES, V5_LB, V5_UB, THETA_NAMES
from isafo.anchor import ANCHOR_LB, ANCHOR_UB, SATURATION_ALARM
from isafo.algos import METHODS, DEFAULT_POP
from isafo.problem import (UNCERTAIN, REL_SPAN, SEED_SPLIT, SEED_TEST,
                           SEED_BASE, split_design_validation)

SOURCES = [
    "isafo/params.py", "isafo/oustaloup.py", "isafo/_kernel.py",
    "isafo/evaluate.py", "isafo/folqi.py", "isafo/anchor.py",
    "isafo/deb.py", "isafo/algos.py", "isafo/problem.py",
    "scripts/run_campaign.py", "isafo/stats.py",
]

N_SEEDS = 30
BUDGET_STAGE1 = 700      # par.3.3: 600 a 800 appels
BUDGET_STAGE2 = 300      # par.3.3: 300 appels minimum
MIN_FEASIBLE_STAGE1 = 20 # critere de passage; sinon graine declaree STERILE
N_TEST_DRAWS = 400       # amendement A1: 381 sont necessaires (Wilson exact)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(65536), b""):
            h.update(blk)
    return h.hexdigest()


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    design, validation = split_design_validation(SEED_SPLIT)

    proto = {
        "nom": "campagne-v6-folqi-fopidf",
        "fige_le": datetime.datetime.now(datetime.timezone.utc)
                    .replace(microsecond=0).isoformat(),
        "objet": ("Optimisation des parametres FO-LQI par six metaheuristiques, "
                  "puis conversion vers FO-PIDF sur Buck non ideal d'ordre entier."),

        "invariants": invariants_dict(),

        "evaluateur": {
            "type": "modele COMMUTE avec blocage reel de la diode",
            "justification": (
                "L'espace de travail du modele Simulink de reference documente "
                "que le palier 1 V siege a 3 % de la frontiere CCM/DCM et que le "
                "saut 24 -> 1 V force environ 3.2 ms de conduction discontinue. "
                "Un evaluateur moyenne a conduction continue serait donc "
                "structurellement faux la ou vit la contrainte g4."),
            "fronts_pwm": "decoupage exact du pas d'integration sur chaque front",
            "pas_recherche_s": NUM.h_search,
            "pas_confirmation_s": NUM.h_confirm,
            "ecart_mesure_J_1us_vs_02us": 1.3e-3,
        },

        "decision_v5": {
            "noms": list(V5_NAMES),
            "borne_basse": V5_LB.tolist(),
            "borne_haute": V5_UB.tolist(),
            "note_ordre": "lambda en 7e position, mu en 5e: ne jamais intervertir.",
        },

        "correcteur_theta": {
            "noms": list(THETA_NAMES),
            "c": 0.0,
            "ancrage_borne_basse": ANCHOR_LB.tolist(),
            "ancrage_borne_haute": ANCHOR_UB.tolist(),
            "seuil_alarme_saturation": SATURATION_ALARM,
        },

        "decision_par_3_1": {
            "question": ("log10(kp) doit-il devenir une 8e variable de decision, "
                         "ou la borne d'ancrage doit-elle etre elargie ?"),
            "reponse": "AUCUNE DES DEUX",
            "justification": (
                "La proposition 1 est reproduite ici a 5.4e-16: kp = Kv + Ki*Vn/(In*R) "
                "est CALCULE exactement par la projection, il n'est pas un degre de "
                "liberte residuel. L'audit de saturation sur 4096 tirages Sobol de "
                "toute la boite v5 montre que log10(kp) ne touche jamais ses bornes "
                "(p05 = -0.64, mediane = +0.44, p95 = +0.96), tres loin de lb = -4. "
                "La saturation a 79 % rapportee au par.2 ne peut donc pas provenir "
                "d'une borne trop serree: il faudrait kp < 1e-4, soit quatre ordres "
                "de grandeur sous ce que produit la projection. Elle signale un "
                "defaut en AMONT de l'ecretage. Promouvoir kp en variable de "
                "decision affaiblirait la part issue de la synthese (le cadrage le "
                "dit lui-meme) pour corriger un symptome et non la cause."),
            "ablation_prevue": ("campagne secondaire avec log10(kp) libere en 8e "
                                "composante, pour mesurer ce que la liberte "
                                "supplementaire apporte reellement"),
        },

        "algorithmes": {
            "methodes": list(METHODS),
            "population": DEFAULT_POP,
            "budget_appariement": "graines communes, budget identique en APPELS",
            "comparateur": "Deb (faisable > infaisable; puis marge; puis violation)",
            "classement_faisables": "-marge (par.3.4), et non J",
        },

        "recherche_etagee": {
            "etage_1": {"cas": "nominal seul", "budget": BUDGET_STAGE1,
                        "critere_passage": f">= {MIN_FEASIBLE_STAGE1} points faisables distincts",
                        "sinon": "graine declaree STERILE et rapportee comme telle"},
            "etage_2": {"cas": "jeu de conception", "budget": BUDGET_STAGE2,
                        "amorcage": "points faisables de l'etage 1 de la MEME graine"},
            "etage_3": {"cas": "replay fin", "pas_s": NUM.h_confirm,
                        "regle": "aucun candidat faisable sans replay fin effectif"},
        },

        "incertitude": {
            "parametres": list(UNCERTAIN),
            "amplitudes_relatives": REL_SPAN,
            "note": ("CHOIX DECLARE: le cadrage evoque 17 cas et 16 coins sans "
                     "les specifier. 16 = 2^4. rC est retenu parce qu'il fixe a "
                     "lui seul le pole projete wf = 1/(rC C)."),
            "Vin": "fixe a 48 V, jamais varie",
        },

        "jeux_disjoints": {
            "graine_partage": SEED_SPLIT,
            "conception": design,
            "validation": validation,
            "graine_test_scellee": SEED_TEST,
            "n_tirages_test": N_TEST_DRAWS,
            "regle": ("un jeu consulte pour decider redevient un jeu de conception "
                      "et perd toute valeur probante"),
        },

        "statistiques": {
            "unite": "la graine, jamais un scenario ni une iteration",
            "n_graines": N_SEEDS,
            "graines": [SEED_BASE + r for r in range(N_SEEDS)],
            "appariement": "methodes appariees sur graines communes",
            "test": "Wilcoxon signe apparie (non parametrique) + Friedman omnibus",
            "correction_multiplicite": "Holm",
            "rapporte": ["medianes", "dispersions", "tailles d'effet",
                         "incertitudes", "TOUS les echecs y compris graines steriles"],
            "interdits": ["conclure sur des p-values intermediaires",
                          "ne rapporter que les graines favorables",
                          "ajuster un seuil apres observation"],
        },

        "limites_declarees": [
            "L'identite FO-LQI -> FO-PIDF n'est exacte que sur le plant LTI "
            "nominal, ip = 0, saturation inactive et etats coherents. La "
            "projection du pole 2e5 -> 2e4 rad/s injecte a elle seule jusqu'a "
            "63 % d'erreur en haut de bande: c'est mesure, pas suppose.",
            "Les snubbers presents dans le modele Simulink ne sont pas reproduits "
            "(valeurs indisponibles); ils n'entrent pas dans J ni dans les six "
            "contraintes.",
            "Le modele de diode est la loi affine statique du constructeur: ni "
            "recouvrement inverse, ni capacite non lineaire, ni retroaction "
            "thermique.",
            "Une reussite sur un nombre fini de cas n'est jamais un certificat "
            "sur un domaine continu.",
        ],
    }

    src = {}
    for rel in SOURCES:
        p = os.path.join(root, rel)
        src[rel] = sha256(p) if os.path.exists(p) else None
    proto["empreintes_sources"] = src

    v6 = os.path.join(root, "protocol", "protocol_v6.json")
    proto["version"] = "6.1"
    proto["amende_depuis"] = {
        "fichier": "protocol/protocol_v6.json",
        "empreinte": json.load(open(v6))["empreinte_protocole"] if os.path.exists(v6) else None,
    }
    proto["statut_des_calculs_anterieurs"] = (
        "Avant cet amendement, seuls des calculs EXPLORATOIRES ont ete faits "
        "(audit de saturation, depistages Sobol, sonde du levier b sur les "
        "graines 1 a 6). Aucune graine de campagne (2026092410 a 2026092439) "
        "n'a ete evaluee, et le jeu de test scelle n'a jamais ete touche.")
    proto["amendements"] = [
        {"id": "A1", "objet": "taille du jeu de test : 300 -> 400",
         "raison": ("Le cadrage indique qu'environ 300 tirages sans echec portent "
                    "la borne de Wilson a 99 %. Le calcul exact donne 98,74 % pour "
                    "300/300 ; il en faut 381. On retient 400 (99,05 %).")},
        {"id": "A2", "objet": "ajout du bras v5b : b en 8e variable de decision",
         "raison": ("Sonde exploratoire (recherche conjointe sur tout theta, "
                    "1500 appels x 6 methodes) : b libre donne une marge "
                    "nominale de +0,191 contre +0,016 avec b = 1, soit 12 fois "
                    "plus, apres replay fin. b = 1 n'est PAS infaisable : c'est un "
                    "levier de marge, pas une barriere. Liberer b affaiblit la part "
                    "issue de la synthese ; c'est declare comme extension, le bras "
                    "v5 conforme au cadrage reste le bras principal.")},
        {"id": "A3", "objet": "l'ablation kp (option A du par.3.1) devient le bras v5kp",
         "raison": ("Annoncee dans v6 sans implementation ; implementee a "
                    "l'identique de v5b, sur les bornes d'ancrage du par.2. "
                    "Priorite la plus basse des trois bras.")},
        {"id": "A4", "objet": "mise a jour des empreintes sources",
         "raison": ("problem.py et run_campaign.py modifies pour A2 et A3 ; "
                    "aucun changement du modele, de la mission, du cout, des "
                    "contraintes ni des seuils.")},
    ]
    proto["bras"] = {"v5": "principal, conforme au cadrage",
                     "v5b": "extension declaree (A2)",
                     "v5kp": "ablation declaree (A3)"}

    body = json.dumps(proto, indent=2, ensure_ascii=False, sort_keys=True)
    proto["empreinte_protocole"] = hashlib.sha256(body.encode()).hexdigest()

    out = os.path.join(root, "protocol", "protocol_v6_1.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(proto, f, indent=2, ensure_ascii=False, sort_keys=True)

    print(f"Protocole fige -> {out}")
    print(f"Empreinte protocole : {proto['empreinte_protocole']}")
    print(f"\nGraines : {SEED_BASE}..{SEED_BASE+N_SEEDS-1} ({N_SEEDS})")
    print(f"Conception : {design}")
    print(f"Validation : {validation}")
    print(f"Budgets : etage 1 = {BUDGET_STAGE1}, etage 2 = {BUDGET_STAGE2}")
    print("\nEmpreintes des sources :")
    for k, v in src.items():
        print(f"  {v[:16]}...  {k}")


if __name__ == "__main__":
    main()
