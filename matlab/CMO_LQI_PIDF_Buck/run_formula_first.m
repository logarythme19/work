%% RUN_FORMULA_FIRST  Chaîne analytique pour le convertisseur de cmo_params.m
%
% Équilibres avec pertes, frontière CCM/DCM, Gvd(s), synthèse LQI (CARE),
% réalisation exacte en PIDF 2-DOF et contrôle des identités de la Proposition 2.
% Aucune boîte à outils requise.

clear; clc;
here = fileparts(mfilename('fullpath')); addpath(here);
F = cmo_lib();
P = cmo_params();                 % <- modifier ici

fprintf('Convertisseur : Vin = %g V, L = %g uH (rL = %g), C = %g uF (rC = %g), R = %g ohm\n', ...
        P.Vin, 1e6*P.L, P.rL, 1e6*P.C, P.rC, P.R);
fprintf('Ron = %g, rd = %g, Vf = %g V, fs = %g kHz, Ts = %g us\n\n', P.Ron, P.rd, P.Vf, P.fs/1e3, 1e6*P.Ts);

fprintf('%6s %8s %8s %9s %10s %9s\n', 'Vo', 'IL', 'd0', 'dIpp', 'marge CCM', 'Rcrit/R');
for V = [P.Vlow P.Vmid P.Vhigh]
    [IL, d0] = F.equilibrium(P, V); c = F.ccm_screen(P, V);
    fprintf('%6.2f %8.3f %8.4f %9.3f %9.1f%% %9.3f\n', V, IL, d0, c.dIpp, c.margin_pct, c.Rcrit/P.R);
end

[num, den] = F.Gvd_coeffs(P, P.Vmid);
fprintf('\nGvd(s) à %g V : num = [%s], den = [%s]\n', P.Vmid, num2str(num, 6), num2str(den, 6));
w0 = sqrt(den(3)/den(1)); z0 = den(2)/den(1)/(2*w0);
fprintf('w0 = %.1f rad/s, zeta0 = %.4f, zéro ESR = %.4g rad/s\n', w0, z0, 1/(P.rC*P.C));

for tag = {'xi0 (Bryson)', 'xi_opt'}
    if strcmp(tag{1}, 'xi_opt'), xi = P.xi_opt; else, xi = P.lqi.xi0; end
    D = F.lqi(P, xi);
    fprintf('\n%s = [%s]\n', tag{1}, num2str(xi, 8));
    fprintf('  K = [%s], résidu CARE relatif %.2e\n', num2str(D.K, 10), D.care_res_rel);
    fprintf('  pôles : %s\n', num2str(D.poles.', 6));
    fprintf('  PIDF : Kp = %.6g, Ki = %.6g, Kd = %.6g, N = %.6g, b = %g, c = %g, Kb = %.6g\n', D.theta);
    [cr, rel] = F.realization_audit(P, D.K);
    fprintf('  résidu des coefficients %.1e, écart fréquentiel relatif max %.1e\n', cr, rel);
    fprintf('  identité Ki_x/rC - Kv_x = Kd N : %.2e\n', abs(D.K(1)/P.rC - D.K(2) - D.theta(3)*D.theta(4)));
end
N = 1/(P.rC*P.C);
fprintf('\nN = %.4g rad/s < pi/Ts = %.4g : %d (condition de discrétisation)\n', N, pi/P.Ts, N < pi/P.Ts);
