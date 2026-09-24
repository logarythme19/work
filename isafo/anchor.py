"""Ancrage du FO-PIDF projete et AUDIT DE SATURATION DES BORNES (par.2).

Le par.2 du cadrage impose de traiter toute borne de la chaine de synthese
comme un objet mesure, pas comme un detail d'implementation: pour chaque
composante, on rapporte la fraction d'appels qui touchent chaque borne, et une
saturation superieure a ~20 % est un signal d'alarme.

Bornes declarees
----------------
La borne basse est celle citee verbatim au par.2:

    lb = [-4, -4, -5, 0.4, 0.4, log10(500), 0]

et la borne haute de log10(kp) y est donnee egale a +2. Les autres bornes
hautes ne figurent pas dans le cadrage; elles sont choisies ici NON SERRANTES
par construction et declarees comme telles, de sorte qu'une saturation
observee sur ces composantes soit imputable a la synthese et non a un choix
arbitraire de l'auteur. La seule exception est log10(wf), dont la borne haute
vaut log10(2e4): c'est le haut de la bande de projection du par.1, et sa
saturation est structurelle (voir isafo.folqi).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .folqi import THETA_NAMES, Synthesis
from .params import NUM

ANCHOR_LB = np.array([-4.0, -4.0, -5.0, 0.40, 0.40, np.log10(500.0), 0.0])
ANCHOR_UB = np.array([2.0, 6.0, -1.0, 1.00, 1.00, np.log10(NUM.proj_wf_hi), 1.0])

# Seuil d'alarme du par.2.
SATURATION_ALARM = 0.20


@dataclass
class SaturationAudit:
    """Compteurs de saturation par composante, cumules sur des appels."""

    n_calls: int = 0
    lo_hits: np.ndarray = field(default_factory=lambda: np.zeros(7, dtype=int))
    hi_hits: np.ndarray = field(default_factory=lambda: np.zeros(7, dtype=int))
    n_nonfinite: np.ndarray = field(default_factory=lambda: np.zeros(7, dtype=int))
    n_rejected: int = 0

    def record(self, theta_raw: np.ndarray, tol: float = 1e-12) -> None:
        self.n_calls += 1
        for j in range(7):
            v = theta_raw[j]
            if not np.isfinite(v):
                self.n_nonfinite[j] += 1
                self.lo_hits[j] += 1     # un gain non positif sature par le bas
                continue
            if v <= ANCHOR_LB[j] + tol:
                self.lo_hits[j] += 1
            elif v >= ANCHOR_UB[j] - tol:
                self.hi_hits[j] += 1

    def fractions(self):
        n = max(self.n_calls, 1)
        return self.lo_hits / n, self.hi_hits / n

    def classify(self, theta_samples: np.ndarray, tol: float = 1e-9):
        """Distingue saturation STRUCTURELLE et saturation PATHOLOGIQUE.

        Une saturation n'a pas toujours le meme sens. Si la composante est une
        CONSTANTE de la synthese (variance nulle sur l'echantillon), la borne
        ne tronque aucun degre de liberte: elle rend seulement visible une
        identite. C'est le cas de log10(wf), fixe a 1/(rC C) par la
        proposition 1, et de b, fixe a 1 par la meme proposition.

        Si au contraire la composante VARIE et vient malgre tout buter sur une
        borne, alors la borne tronque effectivement la recherche: c'est le cas
        pathologique decrit au par.2.
        """
        out = {}
        lo, hi = self.fractions()
        for j, name in enumerate(THETA_NAMES):
            col = theta_samples[:, j]
            fin = col[np.isfinite(col)]
            spread = float(np.ptp(fin)) if fin.size else 0.0
            touches = lo[j] + hi[j]
            if touches <= SATURATION_ALARM:
                kind = "aucune"
            elif spread <= tol:
                kind = "structurelle"
            else:
                kind = "pathologique"
            out[name] = {
                "etendue_theta_brut": spread,
                "fraction_saturee": float(touches),
                "nature": kind,
            }
        return out

    def alarms(self) -> List[str]:
        """Composantes dont la saturation depasse le seuil du par.2."""
        lo, hi = self.fractions()
        out = []
        for j, name in enumerate(THETA_NAMES):
            if lo[j] > SATURATION_ALARM:
                out.append(f"{name}: borne BASSE touchee a {lo[j]*100:.1f} %")
            if hi[j] > SATURATION_ALARM:
                out.append(f"{name}: borne HAUTE touchee a {hi[j]*100:.1f} %")
        return out

    def to_dict(self) -> dict:
        lo, hi = self.fractions()
        return {
            "n_calls": self.n_calls,
            "n_rejected": self.n_rejected,
            "seuil_alarme": SATURATION_ALARM,
            "par_composante": {
                name: {
                    "borne_basse": float(ANCHOR_LB[j]),
                    "borne_haute": float(ANCHOR_UB[j]),
                    "fraction_basse": float(lo[j]),
                    "fraction_haute": float(hi[j]),
                    "n_non_fini": int(self.n_nonfinite[j]),
                    "alarme": bool(lo[j] > SATURATION_ALARM
                                   or hi[j] > SATURATION_ALARM),
                }
                for j, name in enumerate(THETA_NAMES)
            },
            "alarmes": self.alarms(),
        }


def project(syn: Synthesis,
            audit: Optional[SaturationAudit] = None,
            clip: bool = True):
    """Ancrage: theta brut -> theta borne, en enregistrant la saturation.

    Renvoie (theta, sature) ou `sature` indique qu'au moins une composante a
    ete ecretee. Un candidat ecrete n'est PAS invalide, mais il doit etre
    tracable: c'est exactement l'ecretage silencieux qui a fait echouer la
    campagne anterieure.
    """
    theta_raw = syn.theta_raw
    if audit is not None:
        audit.record(theta_raw)
        if not syn.ok:
            audit.n_rejected += 1

    if not clip:
        return theta_raw.copy(), False

    finite = np.where(np.isfinite(theta_raw), theta_raw, ANCHOR_LB)
    theta = np.clip(finite, ANCHOR_LB, ANCHOR_UB)
    saturated = bool(np.any(np.abs(theta - finite) > 1e-12)
                     or not np.all(np.isfinite(theta_raw)))
    return theta, saturated
