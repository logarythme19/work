#!/usr/bin/env python3
"""Figures du par.9.6, statiques pour l'article (PDF vectoriel + PNG).

Regles de trace (skill dataviz): palette validee (3 premiers emplacements,
toutes paires, mode clair), UN SEUL axe par graphique (tension et courant sont
deux graphiques, jamais un double axe), traits fins, limites en encre neutre
pointillee, identite jamais portee par la couleur seule (etiquettes directes et
styles de trait). Les six methodes sont identifiees par leur POSITION sur l'axe,
pas par six teintes, qui ne passeraient pas les seuils daltonisme toutes paires.
"""

import argparse, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8983", "#e4e3df", "#fcfcfb"
LIMIT = dict(color=INK2, lw=0.9, ls=(0, (4, 3)))

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "legend.frameon": False, "lines.linewidth": 1.6,
})


def save(fig, name, outdir):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}.pdf/.png")


def fig_saturation(outdir):
    p = "results/audit_saturation.json"
    if not os.path.exists(p):
        return
    rep = json.load(open(p))
    names = list(rep["par_composante"])
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    for i, n in enumerate(names):
        c = rep["par_composante"][n]
        d = rep["distribution_theta_brut"][n]
        lb, ub = c["borne_basse"], c["borne_haute"]
        lo_, hi_ = min(lb, d["min"]), max(ub, d["max"])
        span = (hi_ - lo_) or 1.0
        f = lambda v: (v - lo_) / span
        ax.plot([i, i], [f(lb), f(ub)], color=GRID, lw=9, solid_capstyle="butt")
        ax.plot([i, i], [f(d["p05"]), f(d["p95"])], color=S1, lw=3,
                solid_capstyle="round")
        ax.plot(i, f(d["median"]), "o", ms=6, color=S1, mec=SURF, mew=1.5)
        nat = rep["classification"][n]["nature"]
        if nat != "aucune":
            ax.annotate(f"{nat}\n{rep['classification'][n]['fraction_saturee']*100:.0f} % saturee",
                        (i, f(d["median"])), xytext=(-8, -4),
                        textcoords="offset points", ha="right", va="top",
                        fontsize=7, color=INK2)
    ax.set_xticks(range(len(names)), names, rotation=25, ha="right")
    ax.set_ylabel("position dans [borne basse, borne haute]")
    ax.set_ylim(-0.05, 1.12)
    ax.set_title(f"Audit de saturation — {rep['n_echantillons']} tirages Sobol de la boite v5",
                 pad=10)
    fig.text(0.01, -0.13, "gris : intervalle autorise [borne basse, borne haute]    "
             "bleu : p05–p95 et mediane de theta brut (un point au-dessus du gris = valeur brute hors borne)",
             fontsize=7, color=INK2)
    save(fig, "fig_audit_saturation", outdir)


def fig_trajectory(theta, Kb, label, outdir, name):
    from isafo.evaluate import Controller
    from isafo.trace import trace
    tr = trace(Controller.from_theta(np.asarray(theta), Kb=Kb), stride=5)
    t = tr["t"] * 1e3
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 4.6), sharex=True)
    a1.plot(t, tr["ref"], color=INK2, lw=1.0, ls=(0, (2, 2)), label="consigne")
    a1.plot(t, tr["vo"], color=S1, lw=1.2, label="v_o")
    a1.axhline(47.5, **LIMIT)
    a1.text(t[-1], 47.5, " Vmax 47,5 V", va="bottom", ha="right", fontsize=7, color=INK2)
    a1.set_ylabel("tension (V)")
    a1.set_title(f"Mission 24 → 1 → 46 V — {label}")
    a1.legend(loc="upper left", fontsize=7)
    a2.plot(t, tr["iL"], color=S2, lw=1.0)
    a2.axhline(20.0, **LIMIT)
    a2.text(t[-1], 20.0, " Imax 20 A", va="bottom", ha="right", fontsize=7, color=INK2)
    a2.set_ylabel("courant i_L (A)")
    a2.set_xlabel("temps (ms)")
    for a in (a1, a2):
        for tx in (108, 216):
            a.axvline(tx, color=GRID, lw=0.8)
    save(fig, name, outdir)


def fig_b_lever(outdir):
    p = "results/probe_setpoint_weight.json"
    if not os.path.exists(p):
        return
    rep = json.load(open(p))
    arms = [("b_libre", "b libre (8e degre de liberte)", S1),
            ("b_fixe_1_proposition_1", "b = 1 (proposition 1)", S2)]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    meths = [r["methode"] for r in rep["bras"][arms[0][0]]["lancers"]]
    for k, (key, lab, col) in enumerate(arms):
        m = [r["marge"] for r in rep["bras"][key]["lancers"]]
        x = np.arange(len(m)) + (k - 0.5) * 0.28
        ax.bar(x, m, width=0.26, color=col, label=lab)
    ax.axhline(0, color=INK, lw=0.9)
    ax.text(-0.45, 0.02, "frontiere de faisabilite (m = 0)", va="bottom",
            ha="left", fontsize=7, color=INK2)
    ax.set_xticks(range(len(meths)), meths)
    ax.set_ylabel("meilleure marge nominale m")
    ax.set_title("Levier de ponderation de consigne b — recherche conjointe, 1500 appels")
    ax.legend(loc="lower right", fontsize=7)
    save(fig, "fig_levier_b", outdir)


def fig_campaign(camp_dir, outdir, tag):
    files = [p for p in glob.glob(os.path.join(camp_dir, "*.json"))
             if not os.path.basename(p).startswith("_")]
    if not files:
        return
    from isafo.algos import METHODS
    rows = [json.load(open(p)) for p in files]
    meths = [m for m in METHODS if any(r["methode"] == m for r in rows)]
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for i, m in enumerate(meths):
        v = np.array([r["etage2"]["meilleure_marge"] for r in rows if r["methode"] == m])
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([i, i], [v.min(), v.max()], color=GRID, lw=1.2)
        ax.plot([i, i], [q1, q3], color=S1, lw=6, solid_capstyle="butt")
        ax.plot(i, med, "o", color=SURF, mec=S1, mew=1.6, ms=6)
        jit = (np.random.default_rng(i).random(v.size) - 0.5) * 0.25
        ax.plot(i + 0.32 + jit, v, ".", color=INK2, ms=3, alpha=0.7)
    ax.axhline(0, color=INK, lw=0.9)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_ylim(top=0.3)
    ax.text(-0.45, 0.05, "faisable (m > 0)", fontsize=7, color=INK2, va="bottom")
    ax.set_xticks(range(len(meths)), meths)
    ax.set_ylabel("marge robuste m (symlog : lineaire sur [-1,1])")
    ax.set_title(f"Distribution par graine — bras {tag} (unite = graine)")
    save(fig, f"fig_graines_{tag}", outdir)

    # Frontiere marge / pire J sur les points faisables de l'etage 2
    pts = [(r["etage2"]["meilleure_marge"], r["etage2"]["meilleur_J"])
           for r in rows if r["etage2"]["faisable"]]
    if len(pts) >= 2:
        P = np.array(pts)
        o = np.argsort(-P[:, 0])
        front, bestJ = [], np.inf
        for i in o:
            if P[i, 1] < bestJ:
                front.append(i); bestJ = P[i, 1]
        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        ax.plot(P[:, 0], P[:, 1], "o", color=MUTED, ms=4, label="candidats faisables")
        F = P[front]; F = F[np.argsort(F[:, 0])]
        ax.plot(F[:, 0], F[:, 1], "-o", color=S1, ms=5, label="non domines")
        ax.set_xlabel("marge robuste m (plus grand = plus sur)")
        ax.set_ylabel("pire J (plus petit = meilleur)")
        ax.set_title(f"Frontiere de compromis marge / performance — bras {tag}")
        ax.legend(fontsize=7)
        save(fig, f"fig_frontiere_{tag}", outdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="figures")
    args = ap.parse_args()
    print("Figures:")
    fig_saturation(args.out)
    fig_b_lever(args.out)
    p = "results/probe_setpoint_weight.json"
    if os.path.exists(p):
        rep = json.load(open(p))
        for key, lab, nm in (("b_libre", "meilleur FO-PIDF, b libre", "fig_trajectoire_b_libre"),
                             ("b_fixe_1_proposition_1", "meilleur FO-PIDF, b = 1", "fig_trajectoire_b1")):
            th = rep["bras"][key]["theta_meilleur"]
            fig_trajectory(th[:7], 1.0 / 10 ** th[7], lab, args.out, nm)
    for tag in ("v5", "v5b", "v5kp"):
        fig_campaign(f"results/campaign_{tag}", args.out, tag)


if __name__ == "__main__":
    main()
