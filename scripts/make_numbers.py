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

        # --- campagne confirmatoire, par bras -------------------------------
        import glob
        for arm, tag in (("v5", "Cinq"), ("v5b", "CinqB"), ("v5kp", "CinqKp")):
            rows = [json.load(open(f)) for f in glob.glob(f"results/campaign_{arm}/*.json")
                    if not os.path.basename(f).startswith("_")]
            an = load(f"results/campaign_{arm}/_analyse.json")
            if not rows or not an:
                continue
            M[f"ChiffreN{tag}"] = str(len(rows))
            M[f"ChiffreSteriles{tag}"] = str(sum(r["etage1"]["sterile"] for r in rows))
            M[f"ChiffreNomFais{tag}"] = str(sum(r["etage1"]["n_faisables"] > 0 for r in rows))
            M[f"ChiffreRobFais{tag}"] = str(sum(bool(r["etage2"]["faisable"]) for r in rows))
            M[f"ChiffreMargeEtUn{tag}"] = f"${fmt(max(r['etage1']['meilleure_marge'] for r in rows), 3, lang)}$"
            M[f"ChiffreMargeEtDeux{tag}"] = f"${fmt(max(r['etage2']['meilleure_marge'] for r in rows), 3, lang)}$"
            st = an["statistiques_marge"]
            M[f"ChiffreFriedmanP{tag}"] = sci(st["friedman"]["p"], 1, lang)
            M[f"ChiffreFriedmanChi{tag}"] = f"${fmt(st['friedman']['chi2'], 1, lang)}$"
            M[f"ChiffrePaires{tag}"] = str(sum(v["significatif_holm"] for v in st["paires"].values()))
            med = {m: v["mediane"] for m, v in st["par_methode"].items()}
            best = max(med, key=med.get)
            M[f"ChiffreMeilleure{tag}"] = best.replace("_", "-").replace("AEABC", "pAEABC").replace("IGWO", "pIGWO") if best not in ("PSO", "GA", "ABC") else best
            M[f"ChiffreMedMeilleure{tag}"] = f"${fmt(med[best], 3, lang)}$"
            lo, hi = an["par_methode"][best]["wilson95_faisable"]
            M["ChiffreWilsonHaut"] = pct(hi, 1, lang)
        fl = load("results/plancher_dmin.json")
        if fl:
            M["ChiffrePlancherN"] = str(fl["n_cas_plancher_sous_seuil"])
            M["ChiffrePlancherTot"] = str(len(fl["cas"]))
            M["ChiffrePlancherNominal"] = f"${fmt(fl['cas']['0']['Vmin_plancher_V'], 3, lang)}$~V"
            M["ChiffrePlancherMin"] = f"${fmt(min(v['Vmin_plancher_V'] for v in fl['cas'].values()), 3, lang)}$~V"

        # --- validation croisee Simulink ---------------------------------
        import csv
        pc, ps = "validation_simulink/candidates.csv", "validation_simulink/simulink_results.csv"
        if os.path.exists(pc) and os.path.exists(ps):
            A = {r["nom"]: r for r in csv.DictReader(open(pc))}
            B = {r["nom"]: r for r in csv.DictReader(open(ps))}
            ctl = [n for n in A if n != "feedforward_seul" and n in B]
            dJ = [abs(float(B[n]["J_sl"]) - float(A[n]["J_py"])) / float(A[n]["J_py"]) for n in ctl]
            dg = [max(abs(float(A[n][f"g{i}_py"]) - float(B[n][f"g{i}_sl"])) for i in range(1, 7)) for n in ctl]
            dV = [float(B[n]["Vmin_sl"]) - float(A[n]["Vmin_py"]) for n in ctl if n.startswith("meilleur")]
            M["ChiffreXvalN"] = str(len(ctl))
            M["ChiffreXvalDJmax"] = pct(max(dJ), 2, lang)
            M["ChiffreXvalDgmax"] = f"${fmt(max(dg), 4, lang)}$"
            M["ChiffreXvalDVmin"] = f"${fmt(1e3 * float(np.mean(dV)), 1, lang)}$~mV"
            M["ChiffreXvalVminSl"] = f"${fmt(float(B['meilleur_v5']['Vmin_sl']), 4, lang)}$~V"
            ff = "feedforward_seul"
            M["ChiffreXvalFFdg"] = f"${fmt(abs(float(A[ff]['g4_py']) - float(B[ff]['g4_sl'])), 3, lang)}$"
            M["ChiffreXvalFFde"] = f"${fmt(1e3 * abs(float(A[ff]['e1_py']) - float(B[ff]['e1_sl'])), 1, lang)}$~mV"
            M["ChiffreXvalMargeSl"] = f"${fmt(-max(float(B['sonde_b1'][f'g{i}_sl']) for i in range(1, 7)), 4, lang)}$"

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
