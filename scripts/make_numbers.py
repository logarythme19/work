#!/usr/bin/env python3
"""Genere les macros chiffrees de l'article depuis les resultats JSON.

Aucun nombre de resultat n'est tape a la main dans les sources LaTeX: ils sont
tous produits ici, en deux versions typographiques (virgule decimale en
francais, point en anglais). Si un resultat change, l'article change avec lui.
Une valeur absente produit un marqueur visible [MANQUANT], jamais un nombre
invente.
"""

import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.evaluate import Controller, evaluate
from isafo.folqi import synthesize, folqi_response, fopidf_response
from isafo.params import NUM


def fmt(v, nd, lang):
    s = f"{v:.{nd}f}"
    return s.replace(".", "{,}") if lang == "fr" else s


def sci(v, nd, lang):
    m, e = f"{v:.{nd}e}".split("e")
    m = m.replace(".", "{,}") if lang == "fr" else m
    return f"${m}\\times10^{{{int(e)}}}$"


def pct(v, nd, lang):
    return f"${fmt(100 * v, nd, lang)}~\\%$"


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    proto = load("protocol/protocol_v6_1.json") or load("protocol/protocol_v6.json")
    audit = load("results/audit_saturation.json")
    probe = load("results/probe_setpoint_weight.json")

    # Chiffres recalcules (deterministes, peu couteux)
    syn = synthesize(np.array([0.0, 0.0, 1.0, 1.0, 1.0, -3.0, 1.0]))
    w = np.logspace(0, 6, 600)
    ref = folqi_response(syn, w)
    e_prop1 = np.max(np.abs(ref - fopidf_response(syn.kp, syn.ki, syn.kd, 1, 1,
                                                  syn.wf_raw, w)) / np.abs(ref))
    wb = np.logspace(np.log10(NUM.proj_wf_lo), np.log10(NUM.proj_wf_hi), 400)
    rb = folqi_response(syn, wb)
    e_pole = np.max(np.abs(rb - fopidf_response(syn.kp, syn.ki, syn.kd, 1, 1,
                                                syn.wf, wb)) / np.abs(rb))
    c0 = Controller(kp=1e-9, ki=1e-9, kd=0.0, lam=1.0, mu=1.0, wf=2e4, b=1.0, Kb=1e3)
    j1, j2 = evaluate(c0, h=NUM.h_search).J, evaluate(c0, h=NUM.h_confirm).J
    e_pas = abs(j1 - j2) / abs(j2)

    for lang in ("fr", "en"):
        M = {}
        M["ChiffreErreurPropUn"] = sci(e_prop1, 1, lang)
        M["ChiffreErreurProjPole"] = pct(e_pole, 0, lang)
        M["ChiffreEcartPas"] = pct(e_pas, 2, lang)
        if audit:
            M["ChiffreAuditN"] = str(audit["n_echantillons"])
            d = audit["distribution_theta_brut"]["log10_kp"]
            M["ChiffreKpPcinq"] = f"${fmt(d['p05'], 2, lang)}$"
            M["ChiffreKpMed"] = f"${fmt(d['median'], 2, lang)}$"
            M["ChiffreKpPqs"] = f"${fmt(d['p95'], 2, lang)}$"
            M["ChiffreKpSat"] = pct(audit["par_composante"]["log10_kp"]["fraction_basse"], 1, lang)
        if proto:
            e = proto["recherche_etagee"]
            M["ChiffreBudgetUn"] = str(e["etage_1"]["budget"])
            M["ChiffreBudgetDeux"] = str(e["etage_2"]["budget"])
            M["ChiffreMinFaisables"] = e["etage_1"]["critere_passage"].split()[1]
            M["ChiffreGraineSplit"] = str(proto["jeux_disjoints"]["graine_partage"])
            M["ChiffreGraineTest"] = str(proto["jeux_disjoints"]["graine_test_scellee"])
            M["ChiffreNTest"] = str(proto["jeux_disjoints"]["n_tirages_test"])
            M["ChiffreEmpreinte"] = "\\texttt{" + proto["empreinte_protocole"][:16] + "}"
        if probe:
            for key, pre in (("b_libre", "Libre"), ("b_fixe_1_proposition_1", "Un")):
                a = probe["bras"][key]
                M[f"ChiffreMargeB{pre}"] = f"${fmt(a['meilleure_marge_02us'], 3, lang)}$"
                M[f"ChiffreNFaisB{pre}"] = str(sum(r["marge"] >= 0 for r in a["lancers"]))
                M[f"ChiffreNLancers{pre}"] = str(len(a["lancers"]))

        lines = ["% Genere par scripts/make_numbers.py -- NE PAS EDITER A LA MAIN"]
        for k, v in sorted(M.items()):
            lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")
        # Tout nom attendu mais absent devient un marqueur visible.
        for k in ("ChiffreMargeBLibre", "ChiffreMargeBUn", "ChiffreNFaisBLibre",
                  "ChiffreNFaisBUn", "ChiffreNLancersLibre", "ChiffreNLancersUn"):
            if k not in M:
                lines.append(f"\\providecommand{{\\{k}}}{{\\textbf{{[MANQUANT]}}}}")
        out = f"docs/article/{lang}/chiffres.tex"
        open(out, "w").write("\n".join(lines) + "\n")
        print(f"{out}: {len(M)} macros")


if __name__ == "__main__":
    main()
