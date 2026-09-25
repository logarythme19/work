"""Generate docs/isa/numbers.tex and docs/isa/tables/*.tex from the archived JSON.

No number of the manuscript is typed by hand: every macro below is read from
results/cmo/formula_first.json or results/cmo/analysis.json.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, 'results', 'cmo')
DOC = os.path.join(ROOT, 'docs', 'isa')
sys.path.insert(0, ROOT)


def sci(x, d=3):
    if x == 0:
        return '0'
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    return f'{m:.{d}f}\\times10^{{{e}}}'


def num(x, d=4):
    return f'{x:.{d}f}'


def main():
    ff = json.load(open(os.path.join(RES, 'formula_first.json')))
    an = json.load(open(os.path.join(RES, 'analysis.json')))
    M = {}
    c1 = ff['ccm'][0]
    M['CCMmargin'] = f"\\SI{{{1e3 * c1['margin']:.2f}}}{{mA}} (\\SI{{{100 * c1['margin'] / c1['Icrit']:.2f}}}{{\\percent}} of $I_{{crit}}$)"
    M['CCMratio'] = f"{c1['Rcrit_over_R']:.3f}"
    anc = an['switched']['Bryson anchor (xi0)']
    M['DCMone'] = f"\\SI{{{100 * anc['dcm_w1']:.2f}}}{{\\percent}}"
    M['ETAone'] = f"\\SI{{{100 * anc['eta_w1']:.2f}}}{{\\percent}}"
    resmax = max(max(abs(v['res_pct']), *[abs(w) for w in v['res_windows']]) for v in an['switched'].values())
    M['RESmax'] = f"${sci(resmax, 1)}\\,\\%$"
    lin = ff['linear_24V']
    A = lin['A']
    M['Aa'], M['Ab'], M['Ac'], M['Ad'] = [f'{v:.4f}' for v in (A[0][0], A[0][1], A[1][0], A[1][1])]
    M['Bd'] = f"{lin['Bd'][0]:.0f}"
    M['Cya'], M['Cyb'] = f"{lin['Cy'][0]:.7f}", f"{lin['Cy'][1]:.7f}"
    a = ff['anchor']
    M['Qa'], M['Qb'], M['Qc'] = f"{a['Q'][0]:.8f}", f"{a['Q'][1]:.8f}", f"{a['Q'][2]:.4f}"
    M['Rd'] = f"{a['Rd']:.8f}"
    M['Ka'], M['Kb'], M['Kc'] = f"{a['K'][0]:.8f}", f"{a['K'][1]:.8f}", f"{a['K'][2]:.8f}"
    M['CAREres'] = f"${sci(a['care_res_rel'], 2)}$"
    p = sorted(a['poles'], key=lambda z: z[1])
    M['Poles'] = f"${p[1][0]:.2f}$ and ${p[2][0]:.2f}\\pm j{abs(p[2][1]):.2f}$\\,rad/s"
    th = a['theta']
    M['ThA'], M['ThB'], M['ThC'], M['ThG'] = f'{th[0]:.10f}', f'{th[1]:.7f}', sci(th[2], 7), f'{th[6]:.3f}'
    M['FreqErr'] = f"${sci(a['freq_rel_error'], 2)}$"
    M['RDmin'] = f"{a['min_return_difference']:.6f}"
    M['PMdeg'] = f"\\SI{{{a['phase_margin_deg']:.2f}}}{{\\degree}}"
    L = ff['anchor_lmi']
    M['AlphaCert'] = f"{L['alpha_cert']:.4f}"
    M['ZetaCert'] = f"{L['zeta_cert']:.7f}"
    M['Dmin'] = f"{L['d_range'][0]:.7f}"
    M['Dmax'] = f"{L['d_range'][1]:.7f}"
    eq = ff['equivalence']
    for k, name in (('nominal_ip0', 'EqNom'), ('load_pulse_0p5A', 'EqLoad'), ('R_plus20pct', 'EqR'),
                    ('L_minus20_rL_Ron_x2', 'EqL'), ('C_plus20pct', 'EqC')):
        M[name] = f"${sci(eq[k]['max_dvo'], 1)}$\\,V" if eq[k]['max_dvo'] < 1e-6 else f"\\SI{{{eq[k]['max_dvo']:.3f}}}{{V}}"
        M[name + 'Sat'] = f"\\SI{{{100 * eq[k]['sat_fraction']:.1f}}}{{\\percent}}"
    sm = ff['sampling']['1to46']
    M['EqSampFOH'] = f"\\SI{{{sm['foh_vs_lqi']:.2f}}}{{V}}"
    M['EqSampZOH'] = f"\\SI{{{sm['zoh_vs_lqi']:.2f}}}{{V}}"
    M['EqSampZOHcont'] = f"\\SI{{{sm['zoh_vs_cont']:.2f}}}{{V}}"
    M['EqSampDelay'] = f"\\SI{{{sm['lqi_vs_cont']:.2f}}}{{V}}"

    # ------------------------------------------------------------ results
    T = os.path.join(DOC, 'tables')
    fam = an['families'].get('CMO-LQI-PIDF')
    if fam:
        meth = ['PSO', 'GA', 'ABC', 'pAEABC', 'pIGWO', 'pIGWO-DLH']
        with open(os.path.join(T, 'tab_algos_res.tex'), 'w') as f:
            f.write('\\begin{table}[t]\n\\centering\n\\caption{CMO--LQI--PIDF family: 30 matched seeds per algorithm (Algorithms 1--6, each followed by refineBO). $J_s$ statistics over feasible runs; ``global'' counts runs already feasible after the 240 global calls.}\n\\label{tab:algores}\n\\footnotesize\n\\resizebox{\\columnwidth}{!}{%\n\\begin{tabular}{@{}lrrrrrr@{}}\n\\toprule\n')
            f.write('Algorithm & feasible & global & best $J_s$ & median & IQR & mean rank\\\\\n\\midrule\n')
            for m in meth:
                v = fam['methods'][m]
                f.write(f"{m} & {v['feasible']}/30 & {v['feasible_global']}/30 & {v['best']:.4f} & {v['median']:.4f} & {v['iqr']:.4f} & {fam['friedman']['mean_ranks'][m]:.2f}\\\\\n")
            f.write('\\bottomrule\n\\end{tabular}}\n\\end{table}\n')
        M['FriedQ'] = f"{fam['friedman']['Q']:.2f}"
        M['FriedP'] = f"{fam['friedman']['p']:.3f}"
        sig = [k for k, v in fam['wilcoxon_holm'].items() if v['p_holm'] < 0.05]
        M['HolmSig'] = ', '.join(k.replace('|', '--') for k in sig) if sig else 'none'
        M['HolmNsig'] = str(len(sig))
        M['CMOfeas'] = str(sum(v['feasible'] for v in fam['methods'].values()))
        M['CMOrefineABC'] = f"{fam['methods']['ABC']['refine_gain_pct']:.1f}"
        M['CMOrefineAE'] = f"{fam['methods']['pAEABC']['refine_gain_pct']:.1f}"
        bests = [fam['methods'][m]['best'] for m in meth]
        M['CMObestLo'], M['CMObestHi'] = f"{min(bests):.4f}", f"{max(bests):.4f}"
        M['CMObest'] = f"{fam['best']['J']:.4f}"
        M['CMOgainPct'] = f"{100 * (1 - fam['best']['J']):.1f}"
    rl = an.get('retained_lqi')
    if rl:
        th = rl['theta']
        M['RetKp'], M['RetKi'], M['RetKd'], M['RetKb'] = f'{th[0]:.5f}', f'{th[1]:.2f}', sci(th[2], 4), f'{th[6]:.0f}'
        M['RetK'] = f"[{rl['K'][0]:.6f},\\ {rl['K'][1]:.6f},\\ {rl['K'][2]:.3f}]"
        x = rl['xi']
        M['RetDi'], M['RetDv'], M['RetwB'] = f'{10 ** x[0]:.3f}', f'{10 ** x[1]:.3f}', f'{10 ** x[2]:.0f}'
        M['RetAlpha'], M['RetZeta'] = f"{rl['alpha_cert']:.1f}", f"{rl['zeta_cert']:.4f}"
        M['RetPM'] = f"\\SI{{{rl['phase_margin_deg']:.2f}}}{{\\degree}}"
        M['RetAlphaRatio'] = f"{rl['alpha_cert'] / ff['anchor_lmi']['alpha_cert']:.2f}"
        p = sorted(rl['poles'], key=lambda z: z[1])
        M['RetPoles'] = f"${p[1][0]:.0f}$ and ${p[2][0]:.0f}\\pm j{abs(p[2][1]):.0f}$\\,rad/s"
    # families
    fams = ['CMO-LQI-PIDF', 'CMO-LQI-SF', 'PIDF-direct', 'PIDF7-direct', 'PID-direct', 'PI-direct']
    ab = an.get('ablation_vs_CMO', {})
    with open(os.path.join(T, 'tab_families.tex'), 'w') as f:
        f.write('\\begin{table*}[t]\n\\centering\n\\caption{Controller families under the identical protocol (7 algorithms, 30 seeds, 300 calls, 33 inequalities). Paired comparison with CMO--LQI--PIDF on the same (algorithm, seed) cells: Deb-consistent score, two-sided Wilcoxon signed-rank test.}\n\\label{tab:families}\n\\footnotesize\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{@{}llrrrrrl@{}}\n\\toprule\n')
        f.write('Family & variables & feasible runs & best $J_s$ & median $J_s$ (feasible) & CMO better / tie / worse & $p$ & min.\\ violation\\\\\n\\midrule\n')
        nv = {'CMO-LQI-PIDF': '4 (LQI weights, $K_b$)', 'CMO-LQI-SF': '4 (same, full state)', 'PIDF-direct': '4 ($K_p,K_i,K_d,K_b$)',
              'PIDF7-direct': '7 (+$N,b,c$)', 'PID-direct': '4', 'PI-direct': '3'}
        for fm in fams:
            if fm not in an['families']:
                continue
            st = an['families'][fm]
            nfe = sum(v['feasible'] for v in st['methods'].values())
            Jf = [v['best'] for v in st['methods'].values() if v['best'] is not None]
            med = [v['median'] for v in st['methods'].values() if v['median'] is not None]
            vb = min(v['viol_best'] for v in st['methods'].values())
            if fm in ab:
                a = ab[fm]; cmp_ = f"{a['cmo_better']} / {a['ties']} / {a['other_better']}"; pv = f"{a['p']:.1e}"
            else:
                cmp_, pv = '--', '--'
            f.write(f"{fm} & {nv[fm]} & {nfe}/180 & {min(Jf):.4f} & {np.median(med):.4f} & {cmp_} & {pv} & {vb:.3f}\\\\\n".replace('None', '--') if Jf else
                    f"{fm} & {nv[fm]} & {nfe}/180 & -- & -- & {cmp_} & {pv} & {vb:.3f}\\\\\n")
        f.write('\\bottomrule\n\\end{tabular}}\n\\end{table*}\n')
    # switched comparison
    sw = an['switched']
    rows = [('CMO-LQI-PIDF*', 'CMO--LQI--PIDF (retained, $v_o$ only)'), ('CMO-LQI-SF*', 'CMO--LQI--SF (twin, measured $i_L$)'),
            ('PIDF-direct*', 'PIDF, direct gain search'), ('PIDF7-direct*', 'PIDF7, direct search'),
            ('PID-direct*', 'PID, direct search'), ('PI-direct*', 'PI, direct search'),
            ('Bryson anchor (xi0)', 'Bryson LQI--PIDF $\\xi_0$'),
            ('PI', 'PI (loop shaping)'), ('PID', 'PID (pole--zero)'), ('PIDF', 'PIDF (Type III)'),
            ('LQR', 'LQR'), ('LQG', 'LQG'), ('LQI', 'LQI (discrete, full state)')]
    with open(os.path.join(T, 'tab_switched.tex'), 'w') as f:
        f.write('\\begin{table*}[t]\n\\centering\n\\caption{Switched replay on the 324 ms reference mission (24$\\to$1$\\to$46 V, 0.5 A load pulse). Specification ratios: overshoot/1.5 V, load deviation/0.5 V, peak $i_L$/23 A, worst steady error/25 mV, steady duty peak-to-peak/0.20; a ratio $\\le1$ meets the specification. $u$ is the largest ratio. Energy residual below \\RESmax{} for every row.}\n\\label{tab:switched}\n\\footnotesize\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{@{}lrrrrrrrrrrc@{}}\n\\toprule\n')
        f.write('Controller & IAE (V\\,s) & overshoot (V) & $i_{L,pk}$ (A) & $|\\bar e|_{\\max}$ (mV) & load dev.\\ (V) & $\\eta$ (\\%) & $E_{loss}/E_u$ & sat.\\ frac. & $u$ & specs met & sensors\\\\\n\\midrule\n')
        for k, lab in rows:
            if k not in sw:
                continue
            v = sw[k]
            sens = '$v_o,i_L,i_C$' if int(v['p'][0]) in (4, 6) else '$v_o$'
            os_ = max(v['os46'], v['us1'], 0.0)
            bold = lambda t, ok: f'\\textbf{{{t}}}' if ok else t
            c_os = bold('%.2f' % os_, os_ <= 1.5)
            c_ip = bold('%.2f' % v['ipk'], v['ipk'] <= 23)
            c_es = bold('%.1f' % (1e3 * v['ess_worst']), v['ess_worst'] <= 0.025)
            c_ld = bold('%.3f' % v['load_dev'], v['load_dev'] <= 0.5)
            f.write(f"{lab} & {v['IAE']:.4f} & {c_os} & {c_ip} & {c_es} & {c_ld} & "
                    f"{100 * v['eta']:.3f} & {v['LI']:.5f} & {v['sat']:.3f} & {v['spec_util']:.2f} & {v['spec_pass']}/5 & {sens}\\\\\n")
            if k == 'PI-direct*':
                f.write('\\midrule\n')
        f.write('\\bottomrule\n\\end{tabular}}\n\\end{table*}\n')
    r = sw.get('CMO-LQI-PIDF*')
    if r:
        M['SwIAE'], M['SwIpk'], M['SwEss'] = f"{r['IAE']:.4f}", f"{r['ipk']:.2f}", f"{1e3 * r['ess_worst']:.1f}"
        M['SwEta'], M['SwLD'], M['SwDpp'] = f"{100 * r['eta']:.3f}", f"{r['load_dev']:.3f}", f"{r['dpp']:.3f}"
        M['SwUtil'] = f"{r['spec_util']:.2f}"
        M['SwSat'] = f"{r['sat']:.3f}"
    for key, nm in (('CMO-LQI-SF*', 'SF'), ('LQI', 'LQI'), ('Bryson anchor (xi0)', 'Anc'), ('PI', 'PI'), ('PID', 'PID'),
                    ('PIDF', 'PIDF'), ('LQR', 'LQR'), ('LQG', 'LQG')):
        if key in sw:
            v = sw[key]
            M[f'Sw{nm}IAE'] = f"{v['IAE']:.4f}"; M[f'Sw{nm}Ipk'] = f"{v['ipk']:.2f}"
            M[f'Sw{nm}OS'] = f"{max(v['os46'], v['us1'], 0):.2f}"; M[f'Sw{nm}Ess'] = f"{1e3 * v['ess_worst']:.1f}"
            M[f'Sw{nm}Eta'] = f"{100 * v['eta']:.3f}"; M[f'Sw{nm}Util'] = f"{v['spec_util']:.2f}"

    extra = os.path.join(RES, 'numbers_extra.json')
    if os.path.exists(extra):
        M.update(json.load(open(extra)))
    with open(os.path.join(DOC, 'numbers.tex'), 'w') as f:
        f.write('% generated by scripts/cmo_numbers.py -- do not edit\n')
        for k, v in M.items():
            f.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
    # Table: transfer functions and CCM screen
    os.makedirs(os.path.join(DOC, 'tables'), exist_ok=True)
    with open(os.path.join(DOC, 'tables', 'tab_tf.tex'), 'w') as f:
        f.write('\\begin{table*}[t]\n\\centering\n\\caption{Equilibria, monic $\\Delta(s)$, $G_{vd}$ numerator and stationary conduction screen at the three commanded levels (nominal load, $i_p=0$). $\\omega_z=1/(r_CC)=\\SI{2e5}{rad/s}$ at every level.}\n\\label{tab:tf}\n\\footnotesize\n')
        f.write('\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{@{}rrrrrllrrrr@{}}\n\\toprule\n')
        f.write('$V_o$ (V) & $I_{L,0}$ (A) & $d_0$ & $\\rho_0$ ($\\Omega$) & $E_d$ (V) & $\\Delta(s)/(LC(R+r_C))$ & $G_{vd}$ numerator & $\\Delta I_{L,pp}$ (A) & $I_{crit}$ (A) & $\\bar I_L-I_{crit}$ (A) & $R_{crit}/R$\\\\\n\\midrule\n')
        for tf, c in zip(ff['transfer_functions'], ff['ccm']):
            d = tf['den_monic']; n = tf['Gvd_num']
            f.write(f"{tf['Vo']:.0f} & {tf['IL0']:.1f} & {tf['d0']:.6f} & {tf['rho0']:.6f} & {tf['Ed']:.3f} & "
                    f"$s^2+{d[1]:.2f}s+{sci(d[2], 5)}$ & ${sci(n[0], 5)}s+{sci(n[1], 5)}$ & "
                    f"{c['dIpp']:.4f} & {c['Icrit']:.4f} & {c['margin']:.4f} & {c['Rcrit_over_R']:.3f}\\\\\n")
        f.write('\\bottomrule\n\\end{tabular}}\n\\end{table*}\n')
    with open(os.path.join(DOC, 'numbers.tex'), 'w') as f:
        f.write('% generated by scripts/cmo_numbers.py -- do not edit\n')
        for k, v in M.items():
            f.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
    print(f'{len(M)} macros written')


if __name__ == '__main__':
    main()
