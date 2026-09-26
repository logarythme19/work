function P = cmo_params(varargin)
% CMO_PARAMS  Tous les paramètres du convertisseur, du modèle Simulink et du protocole.
%
%   P = cmo_params()                       valeurs nominales de l'article
%   P = cmo_params('Vin', 36, 'L', 47e-6)  mêmes valeurs, sauf celles indiquées
%
% Modifiez ce fichier (ou passez des paires nom/valeur) pour étudier un autre
% convertisseur. Les grandeurs dérivées (bornes de recherche, point de
% normalisation, demi-ondulation, scénarios, mission) sont recalculées à la fin
% à partir des valeurs choisies : il suffit de changer les données physiques.
%
% Unités SI : V, A, ohm, H, F, s, Hz.

%% 1. Étage de puissance (blocs du modèle .slx)
P.Vin   = 48;       % tension d'entrée constante (bloc "Vin constant 48 V")
P.L     = 100e-6;   % inductance (bloc L)
P.rL    = 0.10;     % résistance série de l'inductance (bloc rl)
P.C     = 100e-6;   % capacité de sortie (bloc C)
P.rC    = 0.05;     % ESR du condensateur (bloc R1)
P.R     = 10;       % charge résistive nominale (bloc R)
P.Ron   = 0.05;     % résistance à l'état passant du MOSFET (IRF540N)
P.Rbody = 0.01;     % diode interne du MOSFET (Rd du bloc IGBT/Diode)
P.Rs    = 1e5;      % snubber résistif aux bornes de l'interrupteur (Cs = inf)
P.rd    = 0.020;    % résistance de la diode de roue libre (STPST10H100SF)
P.Vf    = 0.42;     % seuil de la diode de roue libre

%% 2. Commande numérique et MLI
P.fs    = 50e3;     % fréquence de la porteuse MLI (bloc PWM1)
P.Ts    = 10e-6;    % période d'échantillonnage du correcteur
P.dmin  = 0.02;     % saturation du rapport cyclique
P.dmax  = 0.98;

%% 3. Niveaux de consigne et perturbation de charge
P.Vmid  = 24;       % niveau de synthèse et niveau initial
P.Vlow  = 1;        % niveau bas
P.Vhigh = 46;       % niveau haut (doit rester < Vin)
P.Iload = 0.5;      % amplitude de l'échelon de courant de charge (A)

%% 4. Mission de validation commutée (Simulink et réplique)
P.mission.t_steps = [0.108 0.216];   % instants des échelons Vmid -> Vlow -> Vhigh
P.mission.T       = 0.324;           % durée totale (Tsim)
P.mission.load    = [0.054 0.090];   % début et fin de l'impulsion de charge
P.mission.h       = 0.2e-6;          % pas RK4 de la réplique commutée

%% 5. Cahier des charges (33 inégalités, valeurs absolues)
P.spec.OS_LIM   = 1.5;     % dépassement maximal (V)
P.spec.LOAD_LIM = 0.5;     % écart maximal sous impulsion de charge (V)
P.spec.I_LIM    = 23;      % courant crête admissible dans l'inductance (A)
P.spec.ESS_TOL  = 0.025;   % erreur statique maximale (V)
P.spec.LI_TOL   = 0.005;   % intensité des pertes <= (1 + LI_TOL) x ancre
P.spec.DPP_LIM  = 0.20;    % excursion crête-à-crête du rapport cyclique (écran commuté)
P.spec.FEAS_TOL = 1e-4;    % bande d'acceptation sur la violation cumulée
P.spec.ESS_WIN  = 2e-3;    % fenêtre de calcul de l'erreur statique (s)

%% 6. Synthèse LQI (règle de Bryson) et point d'ancrage xi0
% xi = [log10 Di, log10 Dv, log10 wB, log10 Kb]
P.lqi.Di0  = 15;            % excursion admissible du courant (A)
P.lqi.Dv0  = [];            % excursion de tension ; vide => Vhigh
P.lqi.wB0  = 2*pi*1500;     % bande de l'intégrateur (rad/s)
P.lqi.Kb0  = 2*pi*1500;     % gain d'anti-windup (1/s)
P.lqi.LB   = [0 0 2 2];     % bornes basses de xi ; bornes hautes calculées à la fin
P.lqi.UB12 = [2 2.4];       % bornes hautes de log10 Di et log10 Dv

%% 7. Paramètres optimaux publiés (CMO-LQI-PIDF, J_s = 0.8151)
P.xi_opt = [1.035116185407974 0.5186738260862431 3.8132348280205437 4.968253545724478];

%% 8. Protocole de recherche
P.search.family   = 'CMO-LQI-PIDF';  % ou CMO-LQI-SF, PIDF-direct, PIDF7-direct, PID-direct, PI-direct
P.search.methods  = {'PSO','GA','ABC','pAEABC','pIGWO','pIGWO-DLH'};
P.search.seeds    = 2026091200 + (1:30);
P.search.npop     = 12;
P.search.B_global = 240;             % appels globaux par exécution
P.search.B_local  = 60;              % appels refineBO par exécution
P.search.nsub     = 5;               % sous-pas RK4 par échantillon (modèle moyen)
P.search.h_screen = 1e-6;            % pas RK4 de l'écran commuté
P.search.verbose  = true;
% bornes log10 [Kp Ki Kd] des familles directes : vide => tirées de la boîte LQI
% (4000 tirages). Valeurs de l'article pour le convertisseur nominal :
% [-4 -0.7 -6.7; 0.5 5.2 -4.3]
P.search.gain_box = [];

%% 9. Modèle Simulink de référence
P.slx     = 'buck_fo_lqi_fixed_reference_v6_1_46V_R2022a';
P.ctrlblk = 'FO_LQI_fixed_reference Command/Common sampled controller';
P.pwmblk  = 'FO_LQI_fixed_reference Command/PWM1';

%% Surcharges nom/valeur (champs de premier niveau ou 'spec.I_LIM', etc.)
for k = 1:2:numel(varargin)
    P = setfield_path(P, varargin{k}, varargin{k+1});
end

%% Grandeurs dérivées
if isempty(P.lqi.Dv0), P.lqi.Dv0 = P.Vhigh; end
assert(P.Vhigh < P.Vin, 'cmo_params: Vhigh doit être inférieur à Vin');
P.lqi.xi0 = log10([P.lqi.Di0 P.lqi.Dv0 P.lqi.wB0 P.lqi.Kb0]);
P.lqi.UB  = [P.lqi.UB12 log10(pi/P.Ts) log10(pi/P.Ts)];
P.lqi.DD  = (P.dmax - P.dmin)/2;
P.lqi.Vsyn = P.Vmid;
% demi-ondulation maximale du courant (ajoutée au pic moyen dans g3)
dI = zeros(1,3); V = [P.Vlow P.Vmid P.Vhigh];
for k = 1:3
    IL = V(k)/P.R;
    d = (V(k) + (P.rL + P.rd)*IL + P.Vf) / (P.Vin - (P.Ron - P.rd)*IL + P.Vf);
    dI(k) = V(k)*(1 - d)/(P.L*P.fs);
end
P.spec.HALF_RIPPLE = 0.5*max(dI);
% six scénarios moyennés : T, instants, niveaux, charge [on off amp], remises à zéro
a = P.Vmid; b = P.Vlow; c = P.Vhigh; I = P.Iload;
P.scen = struct( ...
  'name',  {'S0','S1a','S1b','S1c','S2a','S2b'}, ...
  'T',     {0.090, 0.040, 0.040, 0.040, 0.040, 0.040}, ...
  'tr',    {[0.030 0.060], 0.010, 0.010, 0.010, [], []}, ...
  'lv',    {[a b c], [a b], [b c], [c a], a, c}, ...
  'load',  {[0.015 0.025 I], [0 0 0], [0 0 0], [0 0 0], [0.010 0.020 I], [0.010 0.020 I]}, ...
  'resets',{[0 0.015 0.025 0.030 0.060], [0 0.010], [0 0.010], [0 0.010], [0 0.010 0.020], [0 0.010 0.020]});
% écran commuté pendant la recherche : Vmid -> Vlow -> Vhigh -> Vmid en 34 ms
P.screen.tr  = [0.002 0.014 0.024];
P.screen.lv  = [a b c a];
P.screen.T   = 0.034;
P.screen.win = [0.010 0.014; 0.020 0.024; 0.030 0.034];
% mission : niveaux et fenêtres stationnaires (derniers 10 % de chaque palier)
m = P.mission;
P.mission.levels = [a b c];
edges = [0 m.t_steps m.T];
P.mission.win = [edges(2:end) - 0.1*diff(edges); edges(2:end)]';
P.mission.load_amp = I;
end

function S = setfield_path(S, name, val)
parts = strsplit(name, '.');
if numel(parts) == 1
    S.(name) = val;
else
    S.(parts{1}) = setfield_path(S.(parts{1}), strjoin(parts(2:end), '.'), val);
end
end
