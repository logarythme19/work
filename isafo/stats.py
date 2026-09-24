"""Statistiques predeclarees (par.6).

- L'unite statistique est la GRAINE, jamais un scenario ni une iteration.
- Methodes appariees sur graines communes.
- Test omnibus de Friedman, puis comparaisons par paires en Wilcoxon signe
  apparie, corrigees par la methode de HOLM (nommee, comme l'exige le par.6).
- Tailles d'effet et intervalles, pas seulement des p-values.
- Proportions de reussite avec intervalle de WILSON a 95 %, jamais un verdict
  binaire "robuste" (par.5).
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, List

import numpy as np
from scipy import stats


def wilson(k: int, n: int, z: float = 1.959963984540054):
    """Intervalle de Wilson pour une proportion k/n.

    Reference du par.5: 50/50 donne [92.9 %, 100 %], soit jusqu'a 7 % d'echec
    reel encore compatible avec les donnees.
    """
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def n_for_wilson_lower(target: float, z: float = 1.959963984540054) -> int:
    """Plus petit n tel que n/n succes donne une borne basse >= target."""
    n = 1
    while wilson(n, n, z)[0] < target:
        n += 1
    return n


def holm(pvals: Dict[str, float]) -> Dict[str, float]:
    """Correction de Holm (step-down), p-values ajustees monotones."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    adj, running = {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        adj[k] = running
    return adj


def matched_rank_biserial(a: np.ndarray, b: np.ndarray) -> float:
    """Taille d'effet appariee (correlation rang-biseriale), dans [-1, 1]."""
    d = np.asarray(a) - np.asarray(b)
    d = d[d != 0]
    if d.size == 0:
        return 0.0
    r = stats.rankdata(np.abs(d))
    return float((r[d > 0].sum() - r[d < 0].sum()) / r.sum())


def compare_methods(table: Dict[str, np.ndarray], higher_is_better: bool = True,
                    alpha: float = 0.05) -> dict:
    """table[methode] = vecteur par graine (meme ordre de graines partout)."""
    names = list(table)
    M = np.column_stack([np.asarray(table[n], float) for n in names])
    ok = np.all(np.isfinite(M), axis=1)
    M = M[ok]
    out = {"n_graines_completes": int(M.shape[0]),
           "n_graines_ecartees_non_finies": int((~ok).sum()),
           "par_methode": {}}

    for j, n in enumerate(names):
        col = M[:, j]
        q1, q3 = np.percentile(col, [25, 75]) if col.size else (np.nan, np.nan)
        out["par_methode"][n] = {"mediane": float(np.median(col)) if col.size else None,
                                 "q1": float(q1), "q3": float(q3),
                                 "min": float(col.min()) if col.size else None,
                                 "max": float(col.max()) if col.size else None}

    if M.shape[0] >= 3 and len(names) >= 3:
        # Une colonne constante rend Friedman indefini: on le signale au lieu
        # de fabriquer une p-value.
        if np.all(np.ptp(M, axis=0) == 0) and np.ptp(M) == 0:
            out["friedman"] = {"note": "toutes les valeurs identiques: test indefini"}
        else:
            chi2, p = stats.friedmanchisquare(*[M[:, j] for j in range(M.shape[1])])
            k, n = M.shape[1], M.shape[0]
            ff = ((n - 1) * chi2) / (n * (k - 1) - chi2) if n * (k - 1) != chi2 else np.inf
            out["friedman"] = {"chi2": float(chi2), "p": float(p),
                               "iman_davenport_F": float(ff)}

    raw = {}
    effects = {}
    for (i, a), (j, b) in combinations(enumerate(names), 2):
        x, y = M[:, i], M[:, j]
        key = f"{a} vs {b}"
        if np.allclose(x, y):
            raw[key] = 1.0
        else:
            try:
                raw[key] = float(stats.wilcoxon(x, y, zero_method="wilcox").pvalue)
            except ValueError:
                raw[key] = 1.0
        eff = matched_rank_biserial(x, y)
        effects[key] = eff if higher_is_better else -eff

    adj = holm(raw)
    out["paires"] = {k: {"p_brute": raw[k], "p_holm": adj[k],
                         "effet_rang_biserial": effects[k],
                         "significatif_holm": bool(adj[k] < alpha)}
                     for k in raw}
    out["correction"] = "Holm"
    out["test_apparie"] = "Wilcoxon signe"
    return out
