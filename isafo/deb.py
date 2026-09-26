"""Comparateur de faisabilite d'abord (Deb) et enregistrement des candidats.

Regle du par.3.2: faisable > infaisable; entre faisables, objectif plus petit;
entre infaisables, violation plus faible.

L'"objectif" est volontairement parametrable. Le par.3.4 etablit en effet que
le critere "faisable d'abord, puis J minimal" produit des candidats colles a la
frontiere (marge nominale mediane 0.0066, correlation marge/violation robuste
-0.821), c'est-a-dire les pires points de depart possibles pour la robustesse.
Classer les faisables par -marge plutot que par J transforme le meme
comparateur en chasseur de marge, sans changer la regle de Deb elle-meme.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

MARGIN_TIE = 1e-4      # egalite de marge au sens du par.3.4


@dataclass
class Record:
    """Un appel d'evaluateur, avec tout ce qu'il faut pour l'audit."""

    x: np.ndarray
    J: float
    g: np.ndarray
    margin: float
    feasible: bool
    violation: float
    theta: Optional[np.ndarray] = None
    saturated: bool = False
    n_eval: int = 0
    stage: str = ""
    extra: dict = field(default_factory=dict)

    def active_constraint(self) -> int:
        return int(np.argmax(self.g)) + 1


def _key(rec: Record, mode: str) -> float:
    if mode == "J":
        return rec.J
    if mode == "margin":
        return -rec.margin
    if mode == "lex":
        # marge d'abord, J pour departager a marge quasi egale
        return -(np.floor(rec.margin / MARGIN_TIE) * MARGIN_TIE) + 1e-9 * rec.J
    raise ValueError(f"mode de classement inconnu: {mode}")


def deb_better(a: Record, b: Record, mode: str = "margin") -> bool:
    """Vrai si `a` est strictement meilleur que `b` au sens de Deb."""
    if a is None:
        return False
    if b is None:
        return True
    if a.feasible and not b.feasible:
        return True
    if b.feasible and not a.feasible:
        return False
    if a.feasible and b.feasible:
        return _key(a, mode) < _key(b, mode)
    return a.violation < b.violation


def best_of(records, mode: str = "margin") -> Optional[Record]:
    best = None
    for r in records:
        if deb_better(r, best, mode):
            best = r
    return best


def argsort_deb(records, mode: str = "margin") -> np.ndarray:
    """Indices tries du meilleur au pire, au sens de Deb."""
    def sk(i):
        r = records[i]
        return (0 if r.feasible else 1,
                _key(r, mode) if r.feasible else r.violation)
    return np.array(sorted(range(len(records)), key=sk), dtype=int)
