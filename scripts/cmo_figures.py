"""Figures of the ISA Transactions manuscript (docs/isa/figures/*.pdf)."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from cmo.lqi import lqi, XI0  # noqa: E402
from cmo.equiv import sim_cont  # noqa: E402
from cmo.run import plant_vec, mission_sw  # noqa: E402
from cmo.params import PLANT, SAMP, MISSION_STEPS, MISSION_LEVELS, MISSION_T  # noqa: E402
from cmo.algos import METHODS  # noqa: E402

FIG = os.path.join(ROOT, 'docs', 'isa', 'figures')
RES = os.path.join(ROOT, 'results', 'cmo')
# reference categorical palette, fixed order; line styles as secondary encoding
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
LS = ['-', '--', '-.', ':', (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (1, 1)), '-']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d9d8d4'

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
    'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.5, 'lines.linewidth': 1.2,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})


def fig_equivalence():
    D = lqi(XI0)
    pl = plant_vec()
    tr = np.array(MISSION_STEPS); lv = np.array(MISSION_LEVELS)

    def run(law, amp=0.0, plv=pl):
        return sim_cont(law, plv, PLANT.R, D.K, D.theta, D.theta[6], tr, lv, 0.054, 0.090, amp,
                        2.4, 24.0, MISSION_T, 0.05e-6, SAMP.dmin, SAMP.dmax, 20000)
    a0, b0 = run(0), run(1)
    p3 = pl.copy(); p3[1] = 0.8 * PLANT.L; p3[4] = 2 * PLANT.rL; p3[6] = 2 * PLANT.Ron
    a3, b3 = run(0, 0.0, p3), run(1, 0.0, p3)
    a1, b1 = run(0, 0.5), run(1, 0.5)
    p2 = pl.copy(); p2[3] = 1.2 * PLANT.R
    a2, b2 = run(0, 0.0, p2), run(1, 0.0, p2)
    fig, ax = plt.subplots(2, 1, figsize=(3.5, 3.4), sharex=True)
    t = a0[:, 0] * 1e3
    ax[0].plot(t, a0[:, 1], color=PAL[0], label='LQI (measured $i_L$, $v_C$)')
    ax[0].plot(t, b0[:, 1], color=PAL[1], ls='--', label='2-DOF PIDF (measured $v_o$)')
    ax[0].plot(t, np.interp(a0[:, 0], [0, 0.108, 0.108, 0.216, 0.216, 0.324], [24, 24, 1, 1, 46, 46]),
               color=INK2, lw=0.8, ls=':', label='reference')
    ax[0].set_ylabel('$v_o$ (V)')
    ax[0].legend(loc='upper left')
    fl = 1e-16
    for (a, b), lab, c, s in (((a0, b0), 'nominal, $i_p=0$', PAL[0], '-'),
                              ((a3, b3), '$L$ $-20\\%$, $r_L$, $R_{on}$ $\\times 2$', PAL[2], '-.'),
                              ((a1, b1), '0.5 A load pulse', PAL[3], ':'),
                              ((a2, b2), '$R$ $+20\\%$', PAL[4], '--')):
        ax[1].semilogy(t, np.abs(a[:, 1] - b[:, 1]) + fl, color=c, ls=s, label=lab)
    ax[1].set_ylabel('$|v_o^{LQI}-v_o^{PIDF}|$ (V)')
    ax[1].set_xlabel('time (ms)')
    ax[1].set_ylim(1e-17, 1e8)
    ax[1].set_yticks([1e-15, 1e-10, 1e-5, 1])
    ax[1].legend(loc='upper center', ncol=2, columnspacing=1.0, handlelength=2.2)
    fig.savefig(os.path.join(FIG, 'fig_equivalence.pdf'))
    plt.close(fig)


def fig_campaign(an):
    fam = an['families']['CMO-LQI-PIDF']
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.9))
    for i, m in enumerate(METHODS):
        c = np.array(fam['curves'][m], float)
        c[~np.isfinite(c)] = np.nan
        ax[0].plot(np.arange(1, len(c) + 1), c, color=PAL[i], ls=LS[i], label=m)
    ax[0].axvline(240, color=INK2, lw=0.8, ls=':')
    ax[0].text(243, 0.55, 'refineBO', transform=ax[0].get_xaxis_transform(), va='center', fontsize=7, color=INK2)
    ax[0].set_xlabel('counted evaluator calls')
    ax[0].set_ylabel('median best feasible $J_s$')
    ax[0].legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.22))
    runs = {}
    import glob
    for f in glob.glob(os.path.join(RES, 'campaign', '*', '*.json')):
        r = json.load(open(f))
        runs.setdefault(r['family'], []).append(r)
    fams = ['CMO-LQI-PIDF', 'CMO-LQI-SF', 'PIDF-direct', 'PIDF7-direct', 'PID-direct', 'PI-direct']
    fams = [f for f in fams if f in runs]
    data = [[r['J'] for r in runs[f] if r['feasible']] for f in fams]
    bp = ax[1].boxplot(data, vert=True, widths=0.5, patch_artist=True, showfliers=True,
                       medianprops=dict(color=INK), flierprops=dict(markersize=3, markeredgecolor=INK2))
    for i, b in enumerate(bp['boxes']):
        b.set_facecolor(PAL[i]); b.set_alpha(0.35); b.set_edgecolor(PAL[i])
    ax[1].set_xticks(range(1, len(fams) + 1))
    ax[1].set_xticklabels([f.replace('-direct', '\ndirect').replace('CMO-LQI-', 'CMO-LQI\n') for f in fams])
    for i, f in enumerate(fams):
        n = len(runs[f]); k = sum(r['feasible'] for r in runs[f])
        ax[1].text(i + 1, 1.0, f'{k}/{n}', transform=ax[1].get_xaxis_transform(), ha='center', va='bottom',
                   fontsize=6.5, color=INK2)
    ax[1].set_ylabel('final $J_s$ (feasible runs)')
    ax[1].set_yscale('log')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_campaign.pdf'))
    plt.close(fig)


def fig_switched(an):
    sw = an['switched']
    keys = [('CMO-LQI-PIDF*', 'CMO-LQI-PIDF (retained)'), ('PI', 'PI'), ('PID', 'PID'), ('PIDF', 'PIDF'),
            ('LQR', 'LQR'), ('LQG', 'LQG'), ('LQI', 'LQI')]
    traces = {k: mission_sw(np.array(sw[k]['p'])) for k, _ in keys if k in sw}
    fig = plt.figure(figsize=(7.1, 4.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1])
    a0 = fig.add_subplot(gs[0, :])
    zooms = [fig.add_subplot(gs[1, i]) for i in range(3)]
    for i, (k, lab) in enumerate(keys):
        if k not in traces:
            continue
        s = traces[k]
        t = s['t'] * 1e3
        lw = 1.6 if i == 0 else 1.0
        z = 10 if i == 0 else 3
        a0.plot(t, s['vo'], color=PAL[i], ls=LS[i], lw=lw, label=lab, zorder=z)
        for ax, (x0, x1) in zip(zooms, ((107, 113), (215, 219), (213, 223))):
            m = (t >= x0) & (t <= x1)
            y = s['vo'][m] if ax is not zooms[2] else s['iL'][m]
            ax.plot(t[m], y, color=PAL[i], ls=LS[i], lw=lw, zorder=z)
    s = traces[keys[0][0]]
    a0.plot(s['t'] * 1e3, s['r'], color=INK2, lw=0.8, ls=':', label='reference')
    a0.set_ylabel('$v_o$ (V)'); a0.set_xlabel('time (ms)')
    a0.legend(ncol=4, loc='upper left')
    zooms[0].set_title('24 $\\to$ 1 V'); zooms[0].set_ylabel('$v_o$ (V)')
    zooms[1].set_title('1 $\\to$ 46 V'); zooms[1].set_ylabel('$v_o$ (V)')
    zooms[2].set_title('inductor current, 1 $\\to$ 46 V'); zooms[2].set_ylabel('$i_L$ (A)')
    zooms[2].axhline(23, color=PAL[7], lw=0.8, ls='--')
    zooms[2].text(222.8, 23.5, '23 A limit', ha='right', fontsize=6.5, color=INK2)
    for ax in zooms:
        ax.set_xlabel('time (ms)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_switched.pdf'))
    plt.close(fig)


if __name__ == '__main__':
    os.makedirs(FIG, exist_ok=True)
    fig_equivalence()
    an = json.load(open(os.path.join(RES, 'analysis.json')))
    if 'CMO-LQI-PIDF' in an['families']:
        fig_campaign(an)
    fig_switched(an)
    print('figures written to', FIG)
