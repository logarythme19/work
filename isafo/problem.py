"""Definition du probleme: cas nominal, incertitude, jeux disjoints (par.5).

Ensemble d'incertitude — CHOIX DECLARE
--------------------------------------
Le cadrage evoque "17 cas robustes" et "16 coins" sans dire lesquels. 16 = 2^4:
quatre parametres a deux niveaux, plus le nominal, donnent exactement 17. Le
choix des quatre parametres et de leurs amplitudes n'est PAS dicte par le
cadrage; il est fixe ici et publie dans le protocole:

    L  +/- 20 %      C  +/- 20 %      R  +/- 20 %      rC +/- 50 %

L, C et R gouvernent la dynamique dominante et le point de fonctionnement; rC
est retenu en quatrieme parce qu'il fixe a lui seul le pole projete
wf = 1/(rC C), donc l'erreur de projection elle-meme. Vin reste fixe a 48 V
(par.1), et la mission n'est jamais modifiee.

Regle absolue du par.5: un jeu consulte pour decider redevient un jeu de
conception. Le partage conception/validation est donc TIRE AU SORT avec une
graine publiee, jamais choisi apres avoir vu quels cas echouent.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, Optional, Tuple

import numpy as np

from .params import PLANT, Plant
from .folqi import synthesize, V5_LB, V5_UB
from .anchor import project, SaturationAudit
from .evaluate import Controller, evaluate
from .deb import Record

UNCERTAIN = ("L", "C", "R", "rC")
REL_SPAN = {"L": 0.20, "C": 0.20, "R": 0.20, "rC": 0.50}

# Graines publiees AVANT toute evaluation (par.5 et par.6).
SEED_SPLIT = 2026092401      # partage conception / validation
SEED_TEST = 2026092402       # jeu de test scelle
SEED_BASE = 2026092410       # graines de campagne: SEED_BASE + r


def corner_plants() -> List[Plant]:
    """Le nominal plus les 16 coins de l'hypercube d'incertitude."""
    out = [PLANT]
    for bits in range(16):
        kw = {}
        for j, name in enumerate(UNCERTAIN):
            sign = 1.0 if (bits >> j) & 1 else -1.0
            base = getattr(PLANT, name)
            kw[name] = base * (1.0 + sign * REL_SPAN[name])
        out.append(replace(PLANT, **kw))
    return out


def random_plants(n: int, seed: int) -> List[Plant]:
    """Tirages uniformes dans l'hypercube (jeu de test scelle)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        kw = {}
        for name in UNCERTAIN:
            base = getattr(PLANT, name)
            s = REL_SPAN[name]
            kw[name] = base * rng.uniform(1.0 - s, 1.0 + s)
        out.append(replace(PLANT, **kw))
    return out


def split_design_validation(seed: int = SEED_SPLIT,
                            n_design: int = 8) -> Tuple[List[int], List[int]]:
    """Partage TIRE AU SORT des 16 coins. Le nominal (indice 0) est toujours
    en conception; les coins sont repartis par la graine publiee."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(16) + 1
    design = [0] + sorted(perm[:n_design].tolist())
    validation = sorted(perm[n_design:].tolist())
    return design, validation


# Extension DECLAREE "v5b": b devient une 8e variable de decision.
# La proposition 1 fixe b = 1 exactement; liberer b reconnait ouvertement que
# la ponderation de consigne est un degre de liberte du CORRECTEUR et non de la
# synthese, ce qui AFFAIBLIT d'autant la part issue de la synthese. C'est
# l'analogue exact de l'option "kp en 8e variable" du par.3.1, applique au
# levier que les donnees designent effectivement.
V5B_LB = np.concatenate([V5_LB, [0.0]])
V5B_UB = np.concatenate([V5_UB, [1.0]])


@dataclass
class Problem:
    """v5 -> synthese -> ancrage -> mission(s) -> Record.

    `cases` liste les plants evalues. La marge rapportee est le MINIMUM sur les
    cas (min de max, non differentiable: cf. par.3.4), et J est le PIRE J.
    """

    cases: List[Plant]
    h: float = None
    vref_synthese: float = 24.0
    audit: Optional[SaturationAudit] = None
    clip_anchor: bool = True
    free_b: bool = False

    def __call__(self, x: np.ndarray) -> Record:
        x = np.asarray(x, float)
        syn = synthesize(x[:7], vref=self.vref_synthese)
        if not syn.ok:
            return Record(x=np.asarray(x, float), J=1e6,
                          g=np.full(6, 1e6), margin=-1e6,
                          feasible=False, violation=6e6,
                          theta=None, saturated=True,
                          extra={"rejet": syn.reason})

        theta, sat = project(syn, audit=self.audit, clip=self.clip_anchor)
        if self.free_b:
            if x.shape[0] != 8:
                raise ValueError("v5b attend 8 composantes (v5 + b)")
            theta = theta.copy()
            theta[6] = float(np.clip(x[7], 0.0, 1.0))
        ctrl = Controller.from_theta(theta, Kb=syn.Kb)

        worst_J = -np.inf
        min_margin = np.inf
        tot_viol = 0.0
        g_worst = None
        per_case = []
        for P in self.cases:
            r = evaluate(ctrl, h=self.h, plant=P)
            worst_J = max(worst_J, r.J)
            tot_viol += r.violation
            if r.margin < min_margin:
                min_margin = r.margin
                g_worst = r.g.copy()
            per_case.append({"margin": r.margin, "J": r.J,
                             "feasible": r.feasible,
                             "active": r.active_constraint(),
                             "dcm": r.dcm_fraction})

        return Record(x=np.asarray(x, float), J=float(worst_J),
                      g=g_worst, margin=float(min_margin),
                      feasible=bool(min_margin >= 0.0),
                      violation=float(tot_viol),
                      theta=theta, saturated=sat,
                      extra={"par_cas": per_case, "Kb": syn.Kb,
                             "kp": syn.kp, "ki": syn.ki, "kd": syn.kd})


def nominal_problem(**kw) -> Problem:
    return Problem(cases=[PLANT], **kw)


def robust_problem(indices: List[int], **kw) -> Problem:
    allp = corner_plants()
    return Problem(cases=[allp[i] for i in indices], **kw)
