"""Tests de non-regression de la chaine synthese -> projection -> simulation."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from isafo.folqi import synthesize, folqi_response, fopidf_response
from isafo.evaluate import Controller, evaluate
from isafo.params import NUM, PLANT, MISSION


def test_proposition_1_exacte_en_ordre_entier():
    """alpha = lambda = mu = 1 et pole EXACT -> identite a la precision machine.

    C'est le seul cas ou l'identite FO-LQI -> FO-PIDF est exacte. Tout le reste
    de la campagne vit hors de ce domaine.
    """
    syn = synthesize(np.array([0.0, 0.0, 1.0, 1.0, 1.0, -3.0, 1.0]), vref=24.0)
    assert syn.ok
    w = np.logspace(0, 6, 600)
    err = np.max(np.abs(folqi_response(syn, w)
                        - fopidf_response(syn.kp, syn.ki, syn.kd,
                                          1.0, 1.0, syn.wf_raw, w))
                 / np.abs(folqi_response(syn, w)))
    assert err < 1e-12, f"identite rompue: {err:.3e}"


def test_projection_du_pole_est_une_vraie_erreur():
    """La projection 2e5 -> 2e4 rad/s n'est PAS negligeable: on l'exige mesuree."""
    syn = synthesize(np.array([0.0, 0.0, 1.0, 1.0, 1.0, -3.0, 1.0]), vref=24.0)
    w = np.logspace(np.log10(NUM.proj_wf_lo), np.log10(NUM.proj_wf_hi), 400)
    ref = folqi_response(syn, w)
    err = np.max(np.abs(ref - fopidf_response(syn.kp, syn.ki, syn.kd,
                                              1.0, 1.0, syn.wf, w)) / np.abs(ref))
    assert err > 0.1, "une erreur de projection nulle signalerait un bug de bande"


def test_oustaloup_degenere_proprement():
    """gamma = 0 doit donner l'identite exacte (realisation derivee uniforme)."""
    from isafo.oustaloup import FractionalOperator
    op = FractionalOperator(0.0, NUM.Ts)
    w = np.logspace(0, 5, 50)
    assert np.allclose(op.frequency_response(w), 1.0, atol=1e-12)


def test_feedforward_tient_les_trois_paliers():
    """Sans action de retour, le feedforward nominal doit deja cadrer la mission.

    Garde-fou contre une erreur de signe ou d'echelle dans dff, qui rendrait
    toute la campagne ininterpretable.
    """
    c = Controller(kp=1e-9, ki=1e-9, kd=0.0, lam=1.0, mu=1.0,
                   wf=2e4, b=1.0, Kb=1e3)
    r = evaluate(c)
    assert not r.diverged
    assert abs(r.e24) < 1.0, f"palier 24 V derive: e24={r.e24}"
    assert abs(r.e46) < 1.0, f"palier 46 V derive: e46={r.e46}"


def test_convergence_du_pas_dintegration():
    """L'ecart 1 us -> 0.2 us doit rester petit grace au decoupage des fronts."""
    c = Controller(kp=1e-9, ki=1e-9, kd=0.0, lam=1.0, mu=1.0,
                   wf=2e4, b=1.0, Kb=1e3)
    r1 = evaluate(c, h=NUM.h_search)
    r2 = evaluate(c, h=NUM.h_confirm)
    assert abs(r1.J - r2.J) / abs(r2.J) < 0.02
    assert abs(r1.Imax - r2.Imax) / r2.Imax < 1e-3


def test_dcm_est_bien_atteint():
    """Le palier 1 V doit produire de la conduction discontinue.

    L'espace de travail Simulink documente ~3.2 ms de DCM au saut 24 -> 1 V.
    Un evaluateur qui n'en voit jamais serait un evaluateur moyenne deguise.
    """
    c = Controller(kp=1e-9, ki=1e-9, kd=0.0, lam=1.0, mu=1.0,
                   wf=2e4, b=1.0, Kb=1e3)
    r = evaluate(c)
    assert r.dcm_fraction > 0.0, "aucune conduction discontinue detectee"
