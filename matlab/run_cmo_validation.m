function T = run_cmo_validation(slx_path, csv_path)
% RUN_CMO_VALIDATION  Switched replay of frozen controllers on the reference .slx.
%
%   T = run_cmo_validation('buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx', ...
%                          'controllers.csv')
%
% controllers.csv (exported by scripts/cmo_export_controllers.py) has columns
%   name, p1 ... p36   (the parameter vector of cmo/ctrl.py)
% Only the controller is replaced (FunctionName -> cmo_sfun, Parameters -> p);
% plant, PWM, snubber, solver, mission and loggers are those of the .slx.
% The .slx is never modified or saved (Simulink.SimulationInput only).
% Initial state: iL(0) = 2.4 A, vC(0) = 24 V (equilibrium of the first level).
%
% Metrics use the same definitions as cmo/sim.py and scripts/cmo_validate.py:
%   IAE = trapz |r - vo|, peak iL, max vo, overshoot on the 1->46 V step,
%   steady errors = mean(r - vo) on the last 10 % of each level,
%   energy boundary of ENERGY_ACCOUNTING in the model workspace:
%   Pin = Vin iin, Pu = vo^2/R + vo ip,
%   Ploss = rL iL^2 + rC iC^2 + vSW iSW + vD iD + vSW^2/Rs,
%   residual = Ein - Eu - Eloss - dEstored (percent of Ein).
here = fileparts(mfilename('fullpath'));
addpath(here);
[~, mdl] = fileparts(slx_path);
load_system(slx_path);
cleanup = onCleanup(@() close_system(mdl, 0));

L = 100e-6; C = 100e-6; R = 10; rL = 0.1; rC = 0.05; Rs = 1e5; Vin = 48;
t1 = 0.108; t2 = 0.216; Tsim = 0.324;
win = [0.0972 0.108; 0.2052 0.216; 0.3132 0.324];
ctrl = [mdl '/FO_LQI_fixed_reference Command/Common sampled controller'];
cand = readtable(csv_path, 'TextType', 'string');
out = table();
for k = 1:height(cand)
    p = cand{k, 2:37};
    in = Simulink.SimulationInput(mdl);
    in = in.setBlockParameter(ctrl, 'FunctionName', 'cmo_sfun');
    in = in.setBlockParameter(ctrl, 'Parameters', mat2str(p, 17));
    in = in.setBlockParameter([mdl '/L'], 'InitialCurrent', '2.4');
    in = in.setBlockParameter([mdl '/C'], 'InitialVoltage', '24');
    in = in.setModelParameter('StopTime', num2str(Tsim, 17));
    fprintf('[%d/%d] %s ... ', k, height(cand), cand.name(k));
    tic; so = sim(in); el = toc;

    t  = so.V_o_log.Time;  vo = so.V_o_log.Data(:);
    iL = rs(so.i_L_log, t); iC = rs(so.i_C_log, t); iin = rs(so.i_in_log, t);
    ip = rs(so.i_p_log, t);
    sw = so.igbt_m_log; dd = so.diode_m_log;
    iSW = rs2(sw, t, 1); vSW = rs2(sw, t, 2);
    iD = rs2(dd, t, 1); vD = rs2(dd, t, 2);
    r = 24 * (t < t1) + 1 * (t >= t1 & t < t2) + 46 * (t >= t2);
    vC = vo - rC * iC;
    Pin = Vin * iin; Pu = vo.^2 / R + vo .* ip;
    Pl = rL * iL.^2 + rC * iC.^2 + vSW .* iSW + vD .* iD + vSW.^2 / Rs;
    Es = 0.5 * L * iL.^2 + 0.5 * C * vC.^2;
    Ein = trapz(t, Pin); Eu = trapz(t, Pu); El = trapz(t, Pl);
    res = 100 * (Ein - Eu - El - (Es(end) - Es(1))) / Ein;
    m = struct();
    m.IAE = trapz(t, abs(r - vo));
    m.ipk = max(iL); m.vmax = max(vo);
    m.os46 = max(vo(t >= t2)) - 46;
    for w = 1:3
        sel = t >= win(w, 1) & t <= win(w, 2);
        m.(sprintf('ess%d', w)) = trapz(t(sel), r(sel) - vo(sel)) / (t(find(sel, 1, 'last')) - t(find(sel, 1)));
        m.(sprintf('eta%d', w)) = trapz(t(sel), Pu(sel)) / trapz(t(sel), Pin(sel));
        m.(sprintf('dcm%d', w)) = trapz(t(sel), double(iL(sel) < 1e-3)) / (win(w, 2) - win(w, 1));
    end
    dsat = so.d_sat_log.Data(:); dpre = so.d_pre_log.Data(:);
    m.sat = mean(abs(dsat - dpre) > 1e-12);
    m.eta = Eu / Ein; m.LI = El / Eu; m.Eu = Eu; m.Eloss = El; m.res_pct = res;
    fprintf('IAE=%.5f ipk=%.2f os46=%.3f eta=%.5f res=%.1e%% (%.0f s)\n', m.IAE, m.ipk, m.os46, m.eta, res, el);
    row = struct2table(m); row.name = cand.name(k); row.seconds = el;
    out = [out; row]; %#ok<AGROW>
end
writetable(out, fullfile(here, 'cmo_validation_simulink.csv'));
T = out;
end

function y = rs(ts, t)
[tt, i] = unique(ts.Time, 'last'); d = ts.Data(:); d = d(i);
y = interp1(tt, d, t, 'linear', 'extrap');
end

function y = rs2(ts, t, col)
[tt, i] = unique(ts.Time, 'last'); d = squeeze(ts.Data);
if size(d, 1) ~= numel(ts.Time), d = d.'; end
d = d(i, col);
y = interp1(tt, d, t, 'linear', 'extrap');
end
