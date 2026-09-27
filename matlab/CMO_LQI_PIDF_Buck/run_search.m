%% RUN_SEARCH  Protocole de recherche des paramètres optimaux (Algorithmes 1-7)
%
% Pour chaque algorithme global (PSO, GA, ABC, pAEABC, pIGWO, pIGWO-DLH) et
% chaque graine : B_global appels globaux puis B_local appels refineBO, sur la
% famille choisie (par défaut CMO-LQI-PIDF : recherche dans l'espace des poids
% LQI, chaque candidat étant réalisé exactement en PIDF 2-DOF). Un appel =
% 6 scénarios moyennés + 33 inégalités, dont 3 mesurées sur la réplique commutée.
%
% Fonctionne pour n'importe quel convertisseur Buck : modifiez cmo_params.m ou
% la ligne P = ... (Vin, L, rL, C, rC, R, Ron, rd, Vf, fs, Ts, niveaux, cahier
% des charges). L'ancre de normalisation xi0 et les bornes sont recalculées.
%
% Durée indicative : 2 à 4 s par appel sous MATLAB, soit 10 à 20 min par
% exécution de 300 appels. La campagne complète de l'article (6 algorithmes x
% 30 graines) se lance en parallèle avec use_parfor = true (Parallel Computing
% Toolbox) ou se répartit sur plusieurs machines (listes de graines disjointes).

clear; clc;
here = fileparts(mfilename('fullpath')); addpath(here);
F = cmo_lib();

%% Réglages
P = cmo_params();                          % <- par ex. cmo_params('Vin', 36, 'L', 47e-6)
family  = P.search.family;                 % 'CMO-LQI-PIDF', 'CMO-LQI-SF', 'PIDF-direct', ...
methods = {'pIGWO'};                       % P.search.methods pour les six algorithmes
seeds   = P.search.seeds(1:2);             % P.search.seeds pour les 30 graines de l'article
use_parfor = false;
% budget réduit pour un premier essai (article : 240 + 60)
% P.search.B_global = 60; P.search.B_local = 20;

%% Évaluateur (normalisation par l'ancre de Bryson xi0)
[dec, lb, ub, x0] = F.family(P, family);
fprintf('Famille %s, %d variables, ancre x0 = [%s]\n', family, numel(lb), num2str(x0, 5));
tic; E = F.evaluator(P, dec); fprintf('Référence calculée en %.1f s\n', toc);
r0 = E.fn(x0);
fprintf('Ancre : J = %.4f, faisable = %d, violation = %.3g\n', r0.J, r0.feasible, r0.viol);

%% Exécutions
jobs = {};
for m = 1:numel(methods)
    for s = 1:numel(seeds), jobs(end+1, :) = {methods{m}, seeds(s)}; end %#ok<SAGROW>
end
nj = size(jobs, 1); runs = cell(nj, 1);
if use_parfor
    parfor j = 1:nj
        runs{j} = one_run(P, F, E, jobs{j,1}, jobs{j,2}, lb, ub, x0);
    end
else
    for j = 1:nj
        runs{j} = one_run(P, F, E, jobs{j,1}, jobs{j,2}, lb, ub, x0);
    end
end

%% Synthèse
T = table();
for j = 1:nj
    b = runs{j}.best;
    T = [T; table(string(jobs{j,1}), jobs{j,2}, b.J, b.feasible, b.viol, ...
                  runs{j}.n_feasible, b.n_eval, string(b.stage), {b.x}, ...
                  'VariableNames', {'method','seed','J','feasible','viol','n_feasible', ...
                  'n_eval','stage','x'})]; %#ok<AGROW>
end
disp(T(:, 1:8));
best = F.best_of(cellfun(@(r) r.best, runs, 'UniformOutput', false));
fprintf('\nMeilleur : J = %.4f (faisable = %d), x = [%s]\n', best.J, best.feasible, num2str(best.x, 10));
p_best = dec(best.x);
fprintf('Vecteur p du correcteur : [%s]\n', num2str(p_best(1:8), 8));
if strncmp(family, 'CMO-LQI', 7)
    D = F.lqi(P, best.x);
    fprintf('PIDF 2-DOF : Kp = %.6g, Ki = %.6g, Kd = %.6g, N = %.6g, b = %g, c = %g, Kb = %.6g\n', D.theta);
end
save(fullfile(here, 'cmo_search_result.mat'), 'P', 'family', 'jobs', 'runs', 'T', 'best', 'p_best');
fprintf('Résultats : cmo_search_result.mat. Pour les rejouer : run_optimal (xi = best.x).\n');

function out = one_run(P, F, E, method, seed, lb, ub, x0)
t = tic;
hist = F.run_method(P, method, E.fn, lb, ub, x0, seed);
out.best = F.best_of(hist);
out.best_global = F.best_of(hist(1:min(P.search.B_global, numel(hist))));
out.n_feasible = sum(cellfun(@(r) r.feasible, hist));
curve = nan(1, numel(hist)); cur = NaN;
for k = 1:numel(hist)
    if hist{k}.feasible && (isnan(cur) || hist{k}.J < cur), cur = hist{k}.J; end
    curve(k) = cur;
end
out.curve = curve; out.wall = toc(t);
out.X = cell2mat(cellfun(@(r) r.x, hist, 'UniformOutput', false)');
out.Jall = cellfun(@(r) r.J, hist); out.feas = cellfun(@(r) r.feasible, hist);
if P.search.verbose
    fprintf('%-10s graine %d : J = %.4f, faisable = %d, %d/%d faisables, %.0f s\n', method, seed, ...
            out.best.J, out.best.feasible, out.n_feasible, numel(hist), out.wall);
end
end
