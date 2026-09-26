function T = run_crossval(slx_path)
% RUN_CROSSVAL  Validation croisee Python <-> Simulink des candidats figes.
%
%   T = run_crossval('chemin/vers/buck_fo_lqi_fixed_reference_v6_1_46V_R2022a.slx')
%
% Rejoue chaque ligne de candidates.csv sur le modele Simulink de reference,
% en remplacant SEULEMENT le correcteur par fopidf_sfun (meme interface que
% buck_comparator_sfun). Le fichier .slx n'est JAMAIS modifie ni sauvegarde :
% tout passe par Simulink.SimulationInput.
%
% J et les six contraintes sont calcules avec les MEMES formules que
% isafo/evaluate.py : quadrature sur la grille Ts = 10 us, extrema sur les
% sorties brutes du solveur, erreurs moyennees sur les fenetres stationnaires.
% Resultat : simulink_results.csv, puis lancer compare.py.
%
% HYPOTHESE A VERIFIER : les valeurs de snubber Rs, Cs ne figurent pas dans
% le .slx ; on prend Rs = 1e5 ohm, Cs = inf (snubber resistif). Remplacez-les
% par vos valeurs si vous les connaissez.

here = fileparts(mfilename('fullpath'));
addpath(here);
[~, mdl] = fileparts(slx_path);
load_system(slx_path);
cleanup = onCleanup(@() close_system(mdl, 0));   % jamais de sauvegarde

% ---------------- invariants du cadrage (par.1) --------------------------
Vin = 48; L = 100e-6; C = 100e-6; R = 10; rL = 0.1; rC = 0.05;
Ron = 0.05; rd = 0.020; Vf = 0.42;
Ts = 10e-6; Vn = 24; Dn = 0.48; dmin = 0.02; dmax = 0.98;
Tsim = 0.324; t1 = 0.108; t2 = 0.216;
Rs = 1e5; Cs = inf;                              % HYPOTHESE (voir en-tete)

Vref_ts  = timeseries([24; 1; 46], [0; t1; t2]);
Iload_ts = timeseries([0; 0.5; 0], [0; 0.054; 0.090]);

ctrl = [mdl '/FO_LQI_fixed_reference Command/Common sampled controller'];
cand = readtable(fullfile(here, 'candidates.csv'), 'TextType', 'string');

out = table();
for k = 1:height(cand)
    c = cand(k, :);
    P = [c.kp c.ki c.kd c.lambda c.mu c.wf c.b c.c c.Kb ...
         Ts Vn Dn dmin dmax R rL rd Ron Vin];

    in = Simulink.SimulationInput(mdl);
    in = in.setBlockParameter(ctrl, 'FunctionName', 'fopidf_sfun');
    in = in.setBlockParameter(ctrl, 'Parameters', mat2str(P, 17));
    % Plant nominal et conditions initiales du cadrage (24 V, 2,4 A)
    in = in.setBlockParameter([mdl '/L'],  'Inductance', num2str(L, 17), ...
                              [mdl '/L'],  'InitialCurrent', '2.4');
    in = in.setBlockParameter([mdl '/C'],  'Capacitance', num2str(C, 17), ...
                              [mdl '/C'],  'InitialVoltage', '24');
    in = in.setBlockParameter([mdl '/R'],  'Resistance', num2str(R, 17));
    in = in.setBlockParameter([mdl '/R1'], 'Resistance', num2str(rC, 17));
    in = in.setBlockParameter([mdl '/rl'], 'Resistance', num2str(rL, 17));
    vars = struct('Tsim', Tsim, 'Vref_ts', Vref_ts, 'Iload_ts', Iload_ts, ...
                  'Ron', Ron, 'Vin', Vin, 'd_min', dmin, 'd_max', dmax);
    f = fieldnames(vars);
    for i = 1:numel(f)
        in = in.setVariable(f{i}, vars.(f{i}), 'Workspace', mdl);
    end
    for nm = {'rd', 'Vf', 'Rs', 'Cs', 'R'}
        in = in.setVariable(nm{1}, eval(nm{1}));
    end
    in = in.setModelParameter('StopTime', num2str(Tsim, 17));

    fprintf('[%d/%d] %s ... ', k, height(cand), c.nom);
    tic; so = sim(in); el = toc;

    vo = so.V_o_log;  iL = so.i_L_log;  d = so.d_sat_log;
    m = metrics(vo.Time, vo.Data(:), iL.Time, iL.Data(:), d.Time, d.Data(:), ...
                Ts, Tsim, t1, t2);
    fprintf('J=%.5f marge=%+.4f Imax=%.3f Vmin=%.4f (%.0f s)\n', ...
            m.J, -max(m.g), m.Imax, m.Vmin, el);

    row = table(c.nom, m.J, m.g(1), m.g(2), m.g(3), m.g(4), m.g(5), m.g(6), ...
        m.Imax, m.Vmax, m.Vmin, m.e24, m.e1, m.e46, el, 'VariableNames', ...
        {'nom','J_sl','g1_sl','g2_sl','g3_sl','g4_sl','g5_sl','g6_sl', ...
         'Imax_sl','Vmax_sl','Vmin_sl','e24_sl','e1_sl','e46_sl','secondes'});
    out = [out; row]; %#ok<AGROW>
end
writetable(out, fullfile(here, 'simulink_results.csv'));
fprintf('-> %s\n', fullfile(here, 'simulink_results.csv'));
T = out;
end

function m = metrics(tv, v, ti, i, td, d, Ts, T, t1, t2)
% Memes definitions que isafo/evaluate.py et isafo/_kernel.py.
N  = round(T/Ts);
tk = (0:N-1)' * Ts;
[tv, iv] = unique(tv, 'last'); v = v(iv);
[ti, ii] = unique(ti, 'last'); i = i(ii);
[td, id] = unique(td, 'last'); d = d(id);
vk = interp1(tv, v, tk, 'linear', 'extrap');
ik = interp1(ti, i, tk, 'linear', 'extrap');
dk = interp1(td, d, tk, 'previous', 'extrap');
r  = 24*(tk < t1) + 1*(tk >= t1 & tk < t2) + 46*(tk >= t2);

err = abs(r - vk) ./ max(r, 1);
J = trapz(tk, err)/T ...
  + 0.01 * trapz(tk, ik.^2)/(T*20^2) ...
  + 0.001 * mean(abs(diff([dk(1); dk])))/0.96;

win = @(a, b) mean(r(tk >= a & tk <= b) - vk(tk >= a & tk <= b));
m.e24 = win(0.0936, 0.1044);
m.e1  = win(0.1980, 0.2124);
m.e46 = win(0.3060, 0.3204);
m.Imax = max(abs(i));          % extrema sur les sorties BRUTES du solveur
m.Vmax = max(v);
m.Vmin = min(v);
m.J = J;
m.g = [m.Imax/20 - 1, (m.Vmax - 46)/1.5 - 1, abs(m.e24)/0.05 - 1, ...
       abs(m.e1)/0.03 - 1, abs(m.e46)/0.10 - 1, (0.8 - m.Vmin)/0.8];
end
