"""Invariants physiques et numeriques de la campagne FO-LQI -> FO-PIDF.

Toutes les valeurs de ce module proviennent du paragraphe 1 du cadrage
(PROMPT_CAMPAGNE_FOLQI_ISA.md) et sont NON NEGOCIABLES. Elles sont figees ici
une seule fois et referencees partout ailleurs; aucun autre module ne doit
redefinir une de ces constantes.

Le convertisseur est un Buck non ideal PHYSIQUEMENT D'ORDRE ENTIER. La
fractionnalite appartient au correcteur, jamais au plant.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Tuple


@dataclass(frozen=True)
class Plant:
    """Buck non ideal, source fixe 48 V."""

    Vin: float = 48.0      # V, source constante
    L: float = 100e-6      # H
    C: float = 100e-6      # F
    R: float = 10.0        # Ohm, charge nominale
    rL: float = 0.1        # Ohm, ESR inductance
    rC: float = 0.05       # Ohm, ESR condensateur
    Ron: float = 0.05      # Ohm, resistance de conduction de l'interrupteur
    rd: float = 0.020      # Ohm, resistance dynamique de la diode
    Vf: float = 0.42       # V, tension de seuil de la diode


@dataclass(frozen=True)
class Mission:
    """Mission imposee 24 -> 1 -> 46 V."""

    v0: float = 24.0
    v1: float = 1.0
    v2: float = 46.0
    t_step1: float = 0.108     # s
    t_step2: float = 0.216     # s
    t_end: float = 0.324       # s

    # Perturbation de charge additionnelle
    ip_amp: float = 0.5        # A
    ip_t0: float = 0.054       # s
    ip_t1: float = 0.090       # s

    # Conditions initiales, EN TOUTES CIRCONSTANCES (y compris charge variee).
    # Ce ne sont donc PAS un equilibre perturbe.
    iL0: float = 2.4           # A
    vC0: float = 24.0          # V

    # Fenetres stationnaires d'evaluation des erreurs (s)
    win_24: Tuple[float, float] = (0.0936, 0.1044)
    win_1: Tuple[float, float] = (0.1980, 0.2124)
    win_46: Tuple[float, float] = (0.3060, 0.3204)


@dataclass(frozen=True)
class Numerics:
    """Echantillonnage, PWM, quadrature, approximation d'Oustaloup."""

    Ts: float = 10e-6                 # s, periode d'echantillonnage du correcteur
    f_pwm: float = 50e3               # Hz
    delay_samples: int = 1            # retard discret d'un echantillon
    duty_min: float = 0.02
    duty_max: float = 0.98

    # Pas d'integration: 1 us en recherche, 0.2 us en confirmation (par.5).
    h_search: float = 1e-6            # s
    h_confirm: float = 0.2e-6         # s
    h_quad: float = 10e-6             # s, grille de quadrature du cout

    # Oustaloup: demi-ordre 3 -> 2*3+1 = 7 sections, bande 1..1e5 rad/s
    oustaloup_half_order: int = 3
    oustaloup_wb: float = 1.0
    oustaloup_wh: float = 1e5
    # Bande de projection du pole de filtrage
    proj_wf_lo: float = 20.0
    proj_wf_hi: float = 2e4
    uniform_derivative: bool = True   # active y compris pour mu = 1

    # Normalisations (par.1)
    T0: float = 1e-3                  # s, temps de reference de la synthese
    I_norm: float = 5.0               # A
    V_norm: float = 24.0              # V
    D_norm: float = 0.48              # duty
    R_lqr: float = 1.0                # ponderation d'entree LQR, figee


@dataclass(frozen=True)
class CostSpec:
    """Cout J et six contraintes g <= 0 (par.1)."""

    # J = int|r-v|/max(r,1) dt / T  +  w_i * int iL^2 dt / (T * I_ref^2)
    #     + w_d * mean|dd| / dd_ref
    w_current: float = 0.01
    i_ref: float = 20.0
    w_duty: float = 0.001
    dduty_ref: float = 0.96
    r_floor: float = 1.0              # max(r, 1) au denominateur

    # Contraintes
    imax_lim: float = 20.0            # g1 = Imax/20 - 1
    vmax_ref: float = 46.0            # g2 = (Vmax - 46)/1.5 - 1
    vmax_band: float = 1.5
    e24_tol: float = 0.05             # g3
    e1_tol: float = 0.03              # g4
    e46_tol: float = 0.10             # g5
    vmin_lim: float = 0.8             # g6 = (0.8 - Vmin)/0.8

    n_constraints: int = 6


PLANT = Plant()
MISSION = Mission()
NUM = Numerics()
COST = CostSpec()


def invariants_dict() -> dict:
    """Vue serialisable des invariants, pour l'empreinte du protocole."""
    return {
        "plant": asdict(PLANT),
        "mission": asdict(MISSION),
        "numerics": asdict(NUM),
        "cost": asdict(COST),
    }
