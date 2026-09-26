%% RUN_OPTIMAL  Tester les paramètres optimaux (CMO-LQI-PIDF) et les 6 comparateurs
%
% 1. Calcule, pour le convertisseur décrit dans cmo_params.m, le gain LQI des
%    poids optimaux publiés et sa réalisation exacte en PIDF 2-DOF.
% 2. Rejoue ce correcteur et les comparateurs PI, PID, PIDF, LQR, LQG, LQI
%    sur la mission Vmid -> Vlow -> Vhigh :
%      - sur la réplique commutée MATLAB (aucune boîte à outils requise) ;
%      - sur le modèle Simulink de référence (Simulink + Specialized Power Systems).
% 3. Affiche un tableau spécification par spécification et trace la mission.
%
% Pour un autre convertisseur, modifiez cmo_params.m ou la ligne P = ... ci-dessous,
% par exemple :  P = cmo_params('Vin', 36, 'L', 47e-6, 'Vhigh', 34);
% Les poids publiés ont été optimisés pour le convertisseur nominal : pour un
% autre convertisseur, lancez d'abord run_search.m et utilisez son résultat.

clear; clc;
here = fileparts(mfilename('fullpath')); addpath(here);
F = cmo_lib();

%% Paramètres
P = cmo_params();                 % <- modifier ici (paires nom/valeur)
xi = P.xi_opt;                    % <- ou le résultat de run_search : load('cmo_search_result.mat','best'); xi = best.x;
use_replica  = true;              % réplique commutée MATLAB (≈ 1 min par correcteur)
use_simulink = ~isempty(ver('simulink'));   % modèle .slx de référence
open_model   = false;             % true : ouvre le .slx configuré, prêt pour Run

%% Synthèse LQI -> PIDF 2-DOF exact
D = F.lqi(P, xi);
fprintf('\n=== CMO-LQI-PIDF pour Vin = %g V, L = %g uH, C = %g uF, R = %g ohm ===\n', ...
        P.Vin, 1e6*P.L, 1e6*P.C, P.R);
fprintf('xi = [%s]\n', num2str(xi, 10));
fprintf('Gain LQI K = [Ki_x Kv_x Kz] = [%s]\n', num2str(D.K, 10));
fprintf('Résidu relatif de la CARE : %.2e\n', D.care_res_rel);
fprintf('Pôles en boucle fermée (rad/s) : %s\n', num2str(D.poles.', 6));
th = D.theta;
fprintf('PIDF 2-DOF : Kp = %.6g, Ki = %.6g, Kd = %.6g, N = %.6g, b = %g, c = %g, Kb = %.6g\n', th);
[cr, rel] = F.realization_audit(P, D.K);
fprintf('Audit de la réalisation exacte : résidu %.1e, écart fréquentiel max %.1e\n', cr, rel);

%% Correcteurs à rejouer
C = F.fixed_comparators(P);
ctrl = [struct('name', 'CMO-LQI-PIDF', 'p', F.dec_cmo(P, xi)), ...
        struct('name', 'CMO-LQI-SF', 'p', F.dec_cmo_sf(P, xi)), C];
names = {ctrl.name};
writetable(array2table(cell2mat({ctrl.p}'), 'RowNames', names), ...
           fullfile(here, 'cmo_controllers.csv'), 'WriteRowNames', true);
fprintf('Vecteurs p écrits dans cmo_controllers.csv\n');

%% Réplique commutée MATLAB
S = P.spec;
if use_replica
    fprintf('\n--- Réplique commutée (pas %.1g s) ---\n', P.mission.h);
    fprintf('%-13s %8s %8s %8s %9s %9s %9s %7s %8s\n', 'correcteur', 'IAE', 'iL_pk', ...
            'dep.haut', 'ess1(mV)', 'ess2(mV)', 'ess3(mV)', 'eta', 'specs');
    R = cell(1, numel(ctrl));
    for k = 1:numel(ctrl)
        o = F.mission_sw(P, ctrl(k).p); R{k} = o;
        [ess, os] = mission_indices(P, o);
        ok = [os <= S.OS_LIM, o.ipk <= S.I_LIM, abs(ess) <= S.ESS_TOL];
        fprintf('%-13s %8.4f %8.2f %8.3f %9.1f %9.1f %9.1f %7.4f %5d/5\n', names{k}, o.IAE, ...
                o.ipk, os, 1e3*ess, o.eta, sum(ok));
    end
    figure('Name', 'Mission sur la réplique commutée');
    subplot(2,1,1); hold on;
    for k = [1 3:numel(ctrl)], plot(R{k}.t, R{k}.vo); end
    plot(R{1}.t, R{1}.r, 'k--'); ylabel('v_o (V)'); grid on;
    legend([names([1 3:end]) {'consigne'}], 'Location', 'best');
    subplot(2,1,2); hold on;
    for k = [1 3:numel(ctrl)], plot(R{k}.t, R{k}.iL); end
    yline(S.I_LIM, 'r--'); ylabel('i_L (A)'); xlabel('t (s)'); grid on;
end

%% Modèle Simulink de référence
if use_simulink
    fprintf('\n--- Modèle Simulink %s ---\n', P.slx);
    T = table();
    for k = 1:numel(ctrl)
        fprintf('[%d/%d] %s ... ', k, numel(ctrl), names{k}); tic;
        m = F.sim_simulink(P, ctrl(k).p, here);
        fprintf('IAE = %.4f, iL_pk = %.2f A, dép. = %.3f V (%.0f s)\n', m.IAE, m.ipk, m.os_high, toc);
        row = struct2table(m); row.name = string(names{k}); T = [T; row]; %#ok<AGROW>
    end
    writetable(T, fullfile(here, 'cmo_validation_simulink.csv'));
    disp(T(:, {'name', 'IAE', 'ipk', 'os_high', 'ess1', 'ess2', 'ess3', 'eta', 'res_pct'}));
end
if open_model
    F.configure_model(P, ctrl(1).p, here);
end

function [ess, os] = mission_indices(P, o)
% erreurs statiques (fenêtres stationnaires) et dépassement sur l'échelon haut
ms = P.mission; ess = zeros(1, 3);
for w = 1:3
    sel = o.t >= ms.win(w,1) & o.t < ms.win(w,2);
    ess(w) = mean(o.r(sel) - o.vo(sel));
end
os = max(o.vo(o.t >= ms.t_steps(2))) - ms.levels(3);
end
