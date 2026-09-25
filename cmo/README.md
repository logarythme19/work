# CMO-LQI-PIDF — formula-first reproduction (fixed 48 V non-ideal Buck, 1–46 V)

| Step | Module / script | Output |
|---|---|---|
| Fig. 1 → Eqs. (1)–(15): equilibria, G_vd, CCM screen | `cmo/plant.py` | `results/cmo/formula_first.json` |
| Bryson LQI, CARE, exact 2-DOF PIDF map (Prop. 1) | `cmo/lqi.py` | idem |
| Structural equivalence LQI ≡ PIDF (Prop. 2) | `cmo/equiv.py` | idem, `fig_equivalence.pdf` |
| Common-P D-region LMI certificate | `cmo/lmi.py` | idem |
| Fig. 2 switched replica (snubber, body diode, DCM, exact PWM edges, energy boundary) | `cmo/sim.py` (`sim_sw`) | — |
| Averaged evaluator + switched transition screen, 39 inequalities | `cmo/evaluate.py` | — |
| Algorithms 1–7 (PSO, GA, ABC, pAEABC, pIGWO, pIGWO-DLH, refineBO) | `cmo/algos.py` | — |
| Families and fixed-rule comparators (PI, PID, PIDF, LQR, LQG, LQI) | `cmo/designs.py` | — |
| Campaign (6 families × 6 methods × 30 seeds × 300 calls) | `scripts/cmo_campaign.py` | `results/cmo/campaign/` |
| Statistics, freeze, switched replay | `scripts/cmo_analyze.py` | `results/cmo/analysis.json`, `controllers.csv` |
| MATLAB replay on the reference .slx | `matlab/cmo_sfun.m`, `matlab/run_cmo_validation.m` | `cmo_validation_simulink.csv` |

Run everything from the repository root with `OMP_NUM_THREADS=1`.
