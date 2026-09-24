"""Synthese FO-LQI et projection vers FO-PIDF.

Lecture des indices numeriques du cadrage
-----------------------------------------
Deux nombres du par.1 fixent l'interpretation de toute cette etape.

1. "bande de projection 20-2e4 rad/s". Une bande de projection n'a de sens que
   si quelque chose y est projete: c'est le pole de filtrage derive.
2. Le par.4.3 observe que log10(wf) est "sature a sa borne haute" et juge cette
   borne "justifiee". Or la proposition 1 (article anterieur) donne exactement
   N_C = 1/(rC*C) = 2.0e5 rad/s, tres au-dessus de la bande. La projection de
   ce pole sur [20, 2e4] SATURE donc necessairement, par construction et non
   par accident d'optimisation.

Projection (proposition 1, generalisee a l'ordre fractionnaire)
---------------------------------------------------------------
Pour le plant LTI nominal, ip = 0 et saturation inactive, la relation
entree-sortie du retour d'etat complet se reduit algebriquement a

    C_FOLQI(s) = Kv + a/R + C (a - rC Kv) * N s/(s + N) - Kz/s^alpha
    a = Ki * Vn/In,      N = 1/(rC C)

d'ou l'ancrage FO-PIDF

    kp = Kv + a/R,   kd = C (a - rC Kv),   wf = clip(N, 20, 2e4),
    ki = Kz / T0^alpha,   b = 1,   c = 0.

Le facteur T0^alpha n'est pas cosmetique: l'etat de memoire z est une integrale
d'ordre alpha, donc z_physique = T0^alpha * z_normalise. C'est PAR CE FACTEUR
que l'ordre parent alpha agit reellement sur la synthese, sans ancrage
empirique invente.

L'identite n'est exacte que sur D_eq = {plant LTI nominal, ip = 0, etats
coherents, saturation inactive}. PWM, DCM, anti-windup actif, perturbation de
charge et projection du pole sont HORS de ce domaine. L'erreur de projection
est donc mesuree, jamais supposee nulle (par.4.1).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.linalg import solve_continuous_are

from .params import PLANT, NUM


# Ordre des composantes du vecteur de decision v5 (par.3.1). lambda est en
# SEPTIEME position et mu en CINQUIEME: ne jamais les intervertir.
V5_NAMES = ("log_qi", "log_qv", "log_qmemory", "alpha", "mu",
            "log10_tau_aw_s", "lambda")
V5_LB = np.array([-2.0, -2.0, -1.0, 0.65, 0.40, -4.30, 0.40])
V5_UB = np.array([2.0, 2.0, 3.0, 0.98, 1.00, -2.30, 1.00])

# Ordre des composantes du correcteur (par.4.2).
THETA_NAMES = ("log10_kp", "log10_ki", "log10_kd", "lambda", "mu",
               "log10_wf", "b")


@dataclass
class Linearization:
    A: np.ndarray
    B: np.ndarray
    Cy: np.ndarray
    D_eq: float
    IL_eq: float


def linearize(vref: float, plant=None) -> Linearization:
    """Modele petit-signal normalise du Buck non ideal autour de vref.

    Etats normalises [iL/In, vC/Vn], entree normalisee delta_d/Dn, temps
    normalise t/T0. La sortie vo conserve la respiration de l'ESR.
    """
    P = plant if plant is not None else PLANT
    k = P.R / (P.R + P.rC)

    IL = vref / P.R
    D = (P.R + P.rL + P.rd) * vref / (P.R * P.Vin - (P.Ron - P.rd) * vref)

    a11 = -(P.rd + D * (P.Ron - P.rd) + P.rL + P.rC * k) / P.L
    a12 = -(1.0 - P.rC * k / P.R) / P.L
    a21 = k / P.C
    a22 = -k / (P.R * P.C)
    b1 = (P.Vin - (P.Ron - P.rd) * IL + P.Vf) / P.L

    A_phys = np.array([[a11, a12], [a21, a22]])
    B_phys = np.array([[b1], [0.0]])
    Cy_phys = np.array([[P.rC * k, 1.0 - P.rC * k / P.R]])

    Sx = np.diag([NUM.I_norm, NUM.V_norm])
    Sx_inv = np.diag([1.0 / NUM.I_norm, 1.0 / NUM.V_norm])

    A = NUM.T0 * (Sx_inv @ A_phys @ Sx)
    B = NUM.T0 * (Sx_inv @ B_phys) * NUM.D_norm
    Cy = (Cy_phys @ Sx) / NUM.V_norm

    return Linearization(A=A, B=B, Cy=Cy, D_eq=float(D), IL_eq=float(IL))


@dataclass
class Synthesis:
    Kil: float          # gain sur iL normalise
    Kvc: float          # gain sur vC normalise
    Kz: float           # gain sur l'etat de memoire (temps normalise)
    theta_raw: np.ndarray   # ancrage AVANT tout ecretage
    kp: float
    ki: float
    kd: float
    wf_raw: float
    wf: float
    b: float
    Kb: float
    alpha: float
    lam: float
    mu: float
    ok: bool
    reason: str = ""


def synthesize(x_v5, vref: float = 24.0, plant=None) -> Synthesis:
    """Synthese FO-LQI puis projection FO-PIDF, SANS ecretage.

    L'ecretage eventuel appartient a `anchor.project`, afin que l'audit de
    saturation (par.2) puisse comparer l'ancrage brut a l'ancrage borne.
    """
    x = np.asarray(x_v5, dtype=float)
    if x.shape[0] != 7:
        raise ValueError(f"v5 doit avoir 7 composantes, recu {x.shape[0]}")

    log_qi, log_qv, log_qmem, alpha, mu, log_tau, lam = x

    lin = linearize(vref, plant)
    P = plant if plant is not None else PLANT

    # Augmentation integrale (ordre alpha realise en aval).
    Aa = np.zeros((3, 3))
    Aa[:2, :2] = lin.A
    Aa[2, :2] = -lin.Cy
    Ba = np.zeros((3, 1))
    Ba[:2, 0] = lin.B[:, 0]

    Q = np.diag([10.0 ** log_qi, 10.0 ** log_qv, 10.0 ** log_qmem])
    R = np.array([[NUM.R_lqr]])

    try:
        Pric = solve_continuous_are(Aa, Ba, Q, R)
        K = (np.linalg.solve(R, Ba.T @ Pric)).ravel()
    except Exception as exc:  # CARE non resoluble: candidat rejete, pas corrige
        return Synthesis(Kil=np.nan, Kvc=np.nan, Kz=np.nan,
                         theta_raw=np.full(7, np.nan),
                         kp=np.nan, ki=np.nan, kd=np.nan,
                         wf_raw=np.nan, wf=np.nan, b=1.0, Kb=np.nan,
                         alpha=alpha, lam=lam, mu=mu,
                         ok=False, reason=f"CARE: {exc}")

    Kil, Kvc, Kz = float(K[0]), float(K[1]), float(K[2])

    # --- projection (proposition 1 generalisee) -------------------------
    a = Kil * NUM.V_norm / NUM.I_norm
    kp = Kvc + a / P.R
    kd_time = P.C * (a - P.rC * Kvc)          # en secondes physiques
    N_raw = 1.0 / (P.rC * P.C)

    # Le gain derive est exprime dans le meme temps normalise que le reste
    # du correcteur (le noyau travaille en secondes, donc kd reste en s).
    kd = kd_time
    # Proposition 1: Ki = -Kz. Le signe vient de l'augmentation
    # zdot = -Cy x et de la loi u = -K xa; l'oublier inverse l'action
    # integrale et rend le correcteur inutilisable.
    ki = -Kz / (NUM.T0 ** alpha)

    wf = float(np.clip(N_raw, NUM.proj_wf_lo, NUM.proj_wf_hi))
    Kb = 1.0 / (10.0 ** log_tau)

    with np.errstate(divide="ignore", invalid="ignore"):
        theta_raw = np.array([
            np.log10(kp) if kp > 0 else -np.inf,
            np.log10(ki) if ki > 0 else -np.inf,
            np.log10(kd) if kd > 0 else -np.inf,
            lam, mu,
            np.log10(N_raw),
            1.0,
        ])

    ok = bool(kp > 0 and ki > 0 and kd > 0 and np.all(np.isfinite(K)))
    reason = "" if ok else "gain projete non strictement positif"

    return Synthesis(Kil=Kil, Kvc=Kvc, Kz=Kz, theta_raw=theta_raw,
                     kp=kp, ki=ki, kd=kd, wf_raw=N_raw, wf=wf,
                     b=1.0, Kb=Kb, alpha=alpha, lam=lam, mu=mu,
                     ok=ok, reason=reason)


def folqi_response(syn: Synthesis, omega: np.ndarray, plant=None) -> np.ndarray:
    """Reponse frequentielle EXACTE du retour FO-LQI (avant projection).

    Sert a mesurer l'erreur de projection exigee au par.4.1; c'est la
    reference contre laquelle le FO-PIDF reconstruit est compare.
    """
    P = plant if plant is not None else PLANT
    s = 1j * np.asarray(omega, dtype=complex)
    a = syn.Kil * NUM.V_norm / NUM.I_norm
    feedback = a / P.R + (a * P.C * s + syn.Kvc) / (1.0 + P.rC * P.C * s)
    integral = syn.ki / (s ** syn.alpha)
    return feedback + integral


def fopidf_response(kp, ki, kd, lam, mu, wf, omega) -> np.ndarray:
    """Reponse frequentielle du FO-PIDF reconstruit."""
    s = 1j * np.asarray(omega, dtype=complex)
    return kp + ki / (s ** lam) + kd * (s ** mu) * (wf / (s + wf))
