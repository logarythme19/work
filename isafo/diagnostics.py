"""Diagnostic POST HOC d'oscillation basse frequence (cycle limite).

STATUT : diagnostic ajoute APRES observation. Il est RAPPORTE pour chaque
candidat retenu mais n'intervient JAMAIS dans la selection, le classement ou
la faisabilite, conformement au par.8 (pas de changement du probleme apres
observation).

Motivation : la sonde exploratoire a montre un FO-PIDF faisable au sens des six
contraintes qui entretient un cycle limite a 7,7 kHz au palier 24 V (723 mV
crete a crete). Les erreurs de palier sont des MOYENNES de fenetre, qui
annulent l'oscillation, et le terme en |delta d| de J est trop faible pour la
penaliser.

Indicateur : amplitude crete a crete, dans chaque fenetre stationnaire, de vo
moyennee sur une periode MLI glissante. La moyenne glissante supprime
l'ondulation de decoupage et conserve les oscillations plus lentes.
"""

from __future__ import annotations

import numpy as np

from .params import MISSION, NUM
from .trace import trace


def oscillation_index(ctrl, plant=None, h=None) -> dict:
    tr = trace(ctrl, h=h, plant=plant)
    t, vo, d = tr["t"], tr["vo"], tr["duty"]
    # grille reguliere a 0,5 us puis moyenne glissante sur une periode MLI
    dt = 0.5e-6
    tg = np.arange(t[0], t[-1], dt)
    vg = np.interp(tg, t, vo)
    n = int(round((1.0 / NUM.f_pwm) / dt))
    kern = np.ones(n) / n
    vavg = np.convolve(vg, kern, mode="same")
    out = {}
    for lab, (a, b) in (("24", MISSION.win_24), ("1", MISSION.win_1),
                        ("46", MISSION.win_46)):
        m = (tg >= a + 1e-4) & (tg <= b - 1e-4)
        x = vavg[m] - vavg[m].mean()
        F = np.abs(np.fft.rfft(x * np.hanning(x.size)))
        f = np.fft.rfftfreq(x.size, dt)
        F[f < 200] = 0.0
        out[f"A_osc_{lab}_mV"] = float(np.ptp(vavg[m]) * 1e3)
        out[f"f_dom_{lab}_kHz"] = float(f[np.argmax(F)] / 1e3)
    return out
