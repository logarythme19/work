function F = cmo_lib()
% CMO_LIB  Toutes les fonctions du protocole CMO-LQI-PIDF (MATLAB de base).
%
%   F = cmo_lib();   puis   F.lqi(P, xi), F.evaluator(P, dec), F.run_method(...)
%
% Portage fidèle du paquet Python cmo/ de l'article (mêmes équations, mêmes
% scénarios, mêmes 33 inégalités, mêmes 7 algorithmes). Aucune boîte à outils
% n'est requise : CARE et DARE sont résolues ici (Schur ordonné et doublement).
% Seule la validation sur le .slx demande Simulink + Specialized Power Systems.
%
% Sections :
%   A. Modèle non idéal : équilibre, linéarisation, Gvd
%   B. Synthèse LQI (Bryson -> CARE) et réalisation exacte en PIDF 2-DOF
%   C. Correcteurs de comparaison (PI, PID, PIDF, LQR, LQG, LQI) et familles
%   D. Couche de commande échantillonnée (identique à cmo_sfun.m)
%   E. Simulateurs : modèle moyen (recherche) et réplique commutée (validation)
%   F. Évaluateur : 6 scénarios, J_s, 33 inégalités, écran commuté
%   G. Algorithmes 1-7 : PSO, GA, ABC, pAEABC, pIGWO, pIGWO-DLH, refineBO
%   H. Rejeu sur le modèle Simulink de référence
F = struct();
% A
F.plant_vec = @plant_vec; F.equilibrium = @equilibrium; F.linearize = @linearize;
F.Gvd_coeffs = @Gvd_coeffs; F.ccm_screen = @ccm_screen;
% B
F.care = @care_solve; F.dare = @dare_solve; F.bryson = @bryson; F.lqi = @lqi;
F.pidf_map = @pidf_map; F.realization_audit = @realization_audit;
% C
F.base_vector = @base_vector; F.dec_cmo = @dec_cmo; F.dec_cmo_sf = @dec_cmo_sf;
F.pi_rule = @pi_rule; F.pid_rule = @pid_rule; F.lqr_vec = @lqr_vec;
F.lqg_vec = @lqg_vec; F.lqi_vec = @lqi_vec; F.fixed_comparators = @fixed_comparators;
F.family = @family;
% D
F.ctrl_init = @ctrl_init; F.ctrl_step = @ctrl_step; F.feedforward = @feedforward;
% E
F.sim_avg = @sim_avg; F.sim_sw = @sim_sw; F.mission_sw = @mission_sw;
% F
F.scenario_metrics = @scenario_metrics; F.sw_screen = @sw_screen;
F.evaluator = @evaluator; F.score = @score;
% G
F.run_method = @run_method; F.better = @better; F.best_of = @best_of;
% H
F.sim_simulink = @sim_simulink; F.mission_metrics = @mission_metrics;
F.configure_model = @configure_model; F.model_settings = @model_settings;
end

%% ===================================================== A. modèle non idéal
function pl = plant_vec(P)
pl = [P.Vin P.L P.C P.R P.rL P.rC P.Ron P.rd P.Vf P.Rs P.Rbody];
end

function [IL, d0] = equilibrium(P, Vo, Ip)
% Rapport cyclique d'équilibre avec pertes, Eq. (7)
if nargin < 3, Ip = 0; end
IL = Vo/P.R + Ip;
d0 = (Vo + (P.rL + P.rd)*IL + P.Vf) / (P.Vin - (P.Ron - P.rd)*IL + P.Vf);
end

function lin = linearize(P, Vo, Ip)
% Modèle petits signaux, Eqs. (8)-(10)
if nargin < 3, Ip = 0; end
[IL, d0] = equilibrium(P, Vo, Ip);
g = P.R/(P.R + P.rC);
rho0 = P.rL + P.rd + d0*(P.Ron - P.rd);
Ed = P.Vin + P.Vf - (P.Ron - P.rd)*IL;
lin.A  = [-(rho0 + g*P.rC)/P.L, -g/P.L; g/P.C, -g/(P.R*P.C)];
lin.Bd = [Ed/P.L; 0];
lin.Bg = [d0/P.L; 0];
lin.Bp = [g*P.rC/P.L; -g/P.C];
lin.Cy = [g*P.rC, g];
lin.d0 = d0; lin.IL0 = IL; lin.rho0 = rho0; lin.Ed = Ed; lin.gamma = g;
end

function [num, den] = Gvd_coeffs(P, Vo)
% Gvd(s) en forme fermée, Eqs. (11)-(12) (puissances décroissantes de s)
lin = linearize(P, Vo);
num = lin.Ed*P.R*[P.rC*P.C, 1];
den = [P.L*P.C*(P.R + P.rC), P.L + P.C*((P.R + P.rC)*lin.rho0 + P.R*P.rC), P.R + lin.rho0];
end

function s = ccm_screen(P, Vo)
% Frontière CCM/DCM stationnaire, Eq. (15)
[IL, d] = equilibrium(P, Vo);
dI = Vo*(1 - d)/(P.L*P.fs);
s = struct('Vo', Vo, 'dff', d, 'dIpp', dI, 'Icrit', dI/2, 'margin', IL - dI/2, ...
           'margin_pct', 100*(IL - dI/2)/IL, 'Rcrit', 2*P.L*P.fs/(1 - d));
end

%% ============================================ B. LQI et réalisation PIDF
function X = care_solve(A, B, Q, R)
% Équation de Riccati continue par Schur ordonné de l'hamiltonien
n = size(A, 1);
H = [A, -B*(R\B'); -Q, -A'];
[U, S] = schur(H, 'real');
[U, ~] = ordschur(U, S, diag(S) < 0);
X = U(n+1:end, 1:n)/U(1:n, 1:n);
X = real(X + X')/2;
end

function X = dare_solve(A, B, Q, R)
% Équation de Riccati discrète par l'algorithme de doublement (SDA)
n = size(A, 1); I = eye(n);
Ak = A; Gk = B*(R\B'); Hk = Q;
for it = 1:100
    W = I + Gk*Hk;
    A1 = Ak*(W\Ak);
    G1 = Gk + Ak*(W\Gk)*Ak';
    H1 = Hk + Ak'*Hk*(W\Ak);
    done = norm(H1 - Hk, 1) <= 1e-13*max(1, norm(H1, 1));
    Ak = A1; Gk = G1; Hk = H1;
    if done, break; end
end
X = (Hk + Hk')/2;
end

function [Q, Rd] = bryson(P, xi)
Di = 10^xi(1); Dv = 10^xi(2); wB = 10^xi(3);
Q = diag([Di^-2, Dv^-2, (Dv/wB)^-2]);
Rd = P.lqi.DD^-2;
end

function [Aa, Ba, lin] = augmented(P, Vo)
lin = linearize(P, Vo);
Aa = [lin.A, zeros(2,1); -lin.Cy, 0];
Ba = [lin.Bd; 0];
end

function D = lqi(P, xi)
% xi = [log10 Di, log10 Dv, log10 wB, log10 Kb] -> gain LQI -> PIDF 2-DOF exact
[Aa, Ba] = augmented(P, P.lqi.Vsyn);
[Q, Rd] = bryson(P, xi);
X = care_solve(Aa, Ba, Q, Rd);
K = (Rd\(Ba'*X));                     % [Ki_x, Kv_x, Kz]
res = Aa'*X + X*Aa - X*Ba*(Rd\Ba')*X + Q;
D.xi = xi(:)'; D.Q = Q; D.Rd = Rd; D.P = X; D.K = K;
D.care_res = norm(res); D.care_res_rel = norm(res)/norm(Q);
D.poles = eig(Aa - Ba*K);
D.Kb = 10^xi(4);
D.theta = pidf_map(P, K, D.Kb);
end

function th = pidf_map(P, K, Kb)
% Proposition 1 : Kp = Kv + Ki_x/R, Ki = -Kz, Kd = C (Ki_x - rC Kv), N = 1/(rC C), b = 1, c = 0
Kix = K(1); Kvx = K(2); Kz = K(3);
th = [Kvx + Kix/P.R, -Kz, P.C*(Kix - P.rC*Kvx), 1/(P.rC*P.C), 1, 0, Kb];
end

function [coef_res, rel] = realization_audit(P, K)
% résidu des coefficients et écart fréquentiel relatif maximal LQI / PIDF
th = pidf_map(P, K, 1);
Kix = K(1); Kvx = K(2);
lhs = conv(Kix/P.R, [P.rC*P.C 1]) + [Kix*P.C Kvx];
rhs = conv(th(1), [P.rC*P.C 1]) + [th(3)*th(4)*P.rC*P.C 0];
coef_res = max(abs(lhs - rhs));
lin = linearize(P, P.lqi.Vsyn);
w = logspace(-1, 8, 4000); rel = 0;
for k = 1:numel(w)
    s = 1i*w(k);
    Rv = (s*eye(2) - lin.A)\lin.Bd;
    a = (K(1:2)*Rv)/(lin.Cy*Rv) - K(3)/s;
    b = th(1) + th(2)/s + th(3)*th(4)*s/(s + th(4));
    rel = max(rel, abs(a - b)/abs(a));
end
end

%% ================================== C. comparateurs et familles de recherche
function p = base_vector(P)
p = zeros(1, 36);
p(25:33) = [P.R P.rL P.rd P.Ron P.Vin P.Vf P.rC P.C P.Ts];
p(34:35) = [P.dmin P.dmax];
p(36) = P.Vmid;                       % niveau initial lu par cmo_sfun
end

function p = dec_cmo(P, x)
% CMO-LQI-PIDF : excursions de Bryson -> CARE -> PIDF 2-DOF exact (lit v_o seul)
D = lqi(P, x);
p = base_vector(P); p(1) = 3; p(2:8) = D.theta;
end

function p = dec_cmo_sf(P, x)
% CMO-LQI-SF : même gain continu, état mesuré (i_L, v_C)
D = lqi(P, x);
p = base_vector(P); p(1) = 6; p(2) = D.Kb; p(9:11) = D.K;
end

function G = Gjw(P, w)
[n, d] = Gvd_coeffs(P, P.Vmid);
G = polyval(n, 1i*w)./polyval(d, 1i*w);
end

function p = pi_rule(P)
% PI : coupure 0.03 fs, marge de phase 60 deg, retard 1.5 Ts inclus
DELAY = 1.5*P.Ts; pm = 60;
wc = 2*pi*0.03*P.fs;
G = Gjw(P, wc);
phC = (-180 + pm)*pi/180 - angle(G) + wc*DELAY;
wi = -wc*tan(phC);
Kp = 1/(abs(G)*abs(1 - 1i*wi/wc)); Ki = Kp*wi;
p = base_vector(P); p(1) = 1; p(2) = Kp; p(3) = Ki; p(4) = Ki/Kp;
end

function [p, wc] = pid_rule(P, filtered)
% PID : zéros sur la paire LC ; PIDF (Type III) : + pôle du filtre sur le zéro ESR
if nargin < 2, filtered = false; end
DELAY = 1.5*P.Ts; pm = 60;
[n, d] = Gvd_coeffs(P, P.Vmid);
w0 = sqrt(d(3)/d(1)); z0 = (d(2)/d(1))/(2*w0);
Gdc = n(2)/d(3); wz = 1/(P.rC*P.C);
if filtered
    f = @(w) 90 - w*DELAY*180/pi - pm;
    wc = fzero(f, [1, pi/DELAY]);
    Ki = wc/Gdc;
else
    f = @(w) 90 + atan(w/wz)*180/pi - w*DELAY*180/pi - pm;
    wc = fzero(f, [1, pi/DELAY]);
    Ki = wc/(Gdc*abs(1 + 1i*wc/wz));
end
Kp = Ki*2*z0/w0; Kd = Ki/w0^2;
Kb = 1/sqrt((Kp/Ki)*(Kd/Kp));
p = base_vector(P);
if filtered
    p(1) = 3; p(2:8) = [Kp Ki Kd wz 1 0 Kb];
else
    p(1) = 2; p(2:6) = [Kp Ki Kd 1 Kb];
end
end

function [Ad, Bd] = zoh(A, B, Ts)
n = size(A, 1); m = size(B, 2);
E = expm([A B; zeros(m, n + m)]*Ts);
Ad = E(1:n, 1:n); Bd = E(1:n, n+1:end);
end

function [K, Ad, Bd, lin] = lqr_design(P)
lin = linearize(P, P.Vmid);
[Ad, Bd] = zoh(lin.A, lin.Bd, P.Ts);
[Q, R] = bryson(P, P.lqi.xi0);
X = dare_solve(Ad, Bd, Q(1:2,1:2)*P.Ts, R*P.Ts);
K = (R*P.Ts + Bd'*X*Bd)\(Bd'*X*Ad);
end

function p = lqr_vec(P)
K = lqr_design(P);
p = base_vector(P); p(1) = 4; p(9:10) = K;
end

function p = lqg_vec(P)
% LQR + prédicteur de Kalman stationnaire sur v_o seul
q_d = 0.01; r_v = 0.05;
[K, Ad, Bd, lin] = lqr_design(P);
W = q_d^2*(Bd*Bd'); V = r_v^2; Cd = lin.Cy;
S = dare_solve(Ad', Cd', W, V);
Lf = (S*Cd'/(Cd*S*Cd' + V))';
p = base_vector(P); p(1) = 5; p(9:10) = K; p(11:12) = Lf;
p(13:16) = reshape(Ad', 1, []); p(17:18) = Bd'; p(19:20) = Cd;
end

function p = lqi_vec(P, xi)
if nargin < 2, xi = P.lqi.xi0; end
lin = linearize(P, P.Vmid);
[Ad, Bd] = zoh(lin.A, lin.Bd, P.Ts);
Aa = [Ad, zeros(2,1); -P.Ts*lin.Cy, 1]; Ba = [Bd; 0];
[Q, R] = bryson(P, xi);
X = dare_solve(Aa, Ba, Q*P.Ts, R*P.Ts);
K = (R*P.Ts + Ba'*X*Ba)\(Ba'*X*Aa);
p = base_vector(P); p(1) = 6; p(2) = 10^xi(4); p(9:11) = K;
end

function C = fixed_comparators(P)
C = struct('name', {'PI','PID','PIDF','LQR','LQG','LQI'}, 'p', ...
    {pi_rule(P), pid_rule(P, false), pid_rule(P, true), lqr_vec(P), lqg_vec(P), lqi_vec(P)});
end

function [lo, hi] = gain_box(P)
% plage des gains couverte par la boîte LQI (pour les recherches directes)
st = rng_state(); seed_rng(0); Ur = rand(4000, 4); rng_state(st);
LB = P.lqi.LB; UB = P.lqi.UB; th = zeros(4000, 3);
for k = 1:4000
    D = lqi(P, LB + Ur(k, :).*(UB - LB)); th(k, :) = D.theta(1:3);
end
lo = floor(log10(min(th))*10)/10; hi = ceil(log10(max(th))*10)/10;
end

function [lo, hi] = gain_box_cached(P)
if ~isempty(P.search.gain_box), lo = P.search.gain_box(1,:); hi = P.search.gain_box(2,:);
else, [lo, hi] = gain_box(P); end
end

function [dec, lb, ub, x0] = family(P, name)
% décodeur x -> p, bornes et point initial de chaque famille (Table 7)
LB = P.lqi.LB; UB = P.lqi.UB; XI0 = P.lqi.xi0;
switch name
    case 'CMO-LQI-PIDF'
        dec = @(x) dec_cmo(P, x); lb = LB; ub = UB; x0 = XI0;
    case 'CMO-LQI-SF'
        dec = @(x) dec_cmo_sf(P, x); lb = LB; ub = UB; x0 = XI0;
    case 'PIDF-direct'
        [lo, hi] = gain_box_cached(P); lb = [lo LB(4)]; ub = [hi UB(4)];
        th0 = lqi(P, XI0).theta; x0 = [log10(th0(1:3)) XI0(4)];
        N = 1/(P.rC*P.C);
        dec = @(x) setp(base_vector(P), 1:8, [3 10.^x(1:3) N 1 0 10^x(4)]);
    case 'PIDF7-direct'
        [lo, hi] = gain_box_cached(P);
        lb = [lo 3 0 0 LB(4)]; ub = [hi log10(pi/P.Ts) 1 1 UB(4)];
        th0 = lqi(P, XI0).theta; x0 = [log10(th0(1:4)) 1 0 XI0(4)];
        dec = @(x) setp(base_vector(P), 1:8, [3 10.^x(1:4) x(5) x(6) 10^x(7)]);
    case 'PID-direct'
        [lo, hi] = gain_box_cached(P); lb = [lo LB(4)]; ub = [hi UB(4)];
        p0 = pid_rule(P, false); x0 = min(max(log10(p0([2 3 4 6])), lb), ub);
        dec = @(x) setp(base_vector(P), 1:6, [2 10.^x(1:3) 1 10^x(4)]);
    case 'PI-direct'
        [lo, hi] = gain_box_cached(P); lb = [lo(1:2) LB(4)]; ub = [hi(1:2) UB(4)];
        p0 = pi_rule(P); x0 = min(max(log10(p0(2:4)), lb), ub);
        dec = @(x) setp(base_vector(P), 1:4, [1 10.^x(1:3)]);
    otherwise
        error('cmo:family', 'famille inconnue : %s', name);
end
end

function p = setp(p, idx, v)
p(idx) = v;
end

%% ================================== D. couche de commande échantillonnée
% Vecteur p (36) : identique à cmo_sfun.m (voir l'en-tête de ce fichier-là).
% État s (8) : 1 intégrale, 2 filtre / sortie précédente, 3:4 observateur,
%              5 rapport cyclique appliqué, 6 d_pre, 7 entrée précédente du filtre, 8 drapeau.
function d = feedforward(r, p)
R = p(25); rL = p(26); rd = p(27); Ron = p(28); Vin = p(29); Vf = p(30);
IL = r/R;
d = (r + (rL + rd)*IL + Vf)/(Vin - (Ron - rd)*IL + Vf);
end

function s = ctrl_init(p, r0, vo0, iL0, vC0)
s = zeros(8, 1);
d0 = feedforward(r0, p);
s(5) = d0; s(6) = d0;
switch p(1)
    case 3, s(2) = p(7)*r0 - vC0;     % x_F(0) = c r - v_C(0) (Proposition 2)
    case 2, s(2) = -vo0;
    case 5, s(3) = iL0; s(4) = vC0;
end
end

function [s, dnow, ds, dt] = ctrl_step(p, s, r, vo, iL, iC)
code = p(1); Ts = p(33); rC = p(31); R = p(25);
dmin = p(34); dmax = p(35);
dff = feedforward(r, p);
e = r - vo;
dnow = s(5);
u = 0;
switch code
    case 1
        u = p(2)*e + s(1);
    case 2
        yd = -vo;
        u = p(2)*(p(5)*r - vo) + s(1) + p(4)*(yd - s(2))/Ts;
        s(2) = yd;
    case 3
        Kd = p(4); N = p(5); b = p(6); c = p(7);
        w = c*r - vo;
        if p(9) == 1                  % filtre ZOH historique (audit seulement)
            u = p(2)*(b*r - vo) + s(1) + Kd*N*(w - s(2));
            phi = exp(-N*Ts);
            s(2) = phi*s(2) + (1 - phi)*w;
        else                          % maintien triangulaire (FOH)
            if s(8) == 0
                s(8) = 1;
            else
                phi = exp(-N*Ts); gam = 1 - (1 - phi)/(N*Ts);
                s(2) = phi*s(2) + (1 - phi)*s(7) + gam*(w - s(7));
            end
            s(7) = w;
            u = p(2)*(b*r - vo) + s(1) + Kd*N*(w - s(2));
        end
    case {4, 6}
        vC = vo - rC*iC;
        u = -p(9)*(iL - r/R) - p(10)*(vC - r);
        if code == 6, u = u - p(11)*s(1); end
    case 5
        xh0 = s(3); xh1 = s(4); xe0 = r/R; xe1 = r;
        yh = r + p(19)*(xh0 - xe0) + p(20)*(xh1 - xe1);
        nu = vo - yh;
        u = -p(9)*(xh0 - xe0) - p(10)*(xh1 - xe1);
        a0 = xh0 + p(11)*nu - xe0; a1 = xh1 + p(12)*nu - xe1;
        s(3) = xe0 + p(13)*a0 + p(14)*a1 + p(17)*(dnow - dff);
        s(4) = xe1 + p(15)*a0 + p(16)*a1 + p(18)*(dnow - dff);
end
dt = dff + u;
ds = min(max(dt, dmin), dmax);
switch code
    case 1, s(1) = s(1) + Ts*(p(3)*e + p(4)*(ds - dt));
    case 2, s(1) = s(1) + Ts*(p(3)*e + p(6)*(ds - dt));
    case 3, s(1) = s(1) + Ts*(p(3)*e + p(8)*(ds - dt));
    case 6, s(1) = s(1) + Ts*(e + p(2)/abs(p(11))*(ds - dt));
end
s(5) = ds; s(6) = dt;
end

%% ================================================= E. simulateurs
function r = ref(t, tr, lv)
r = lv(1);
for k = 1:numel(tr)
    if t >= tr(k), r = lv(k+1); end
end
end

function [ic, vo] = cap(iL, vC, ip, R, rC)
ic = (iL - vC/R - ip)/(1 + rC/R);
vo = vC + rC*ic;
end

function [diL, dvC, pin, pu, pls] = avg_f(iL, vC, d, ip, pl)
Vin = pl(1); L = pl(2); C = pl(3); R = pl(4); rL = pl(5); rC = pl(6);
Ron = pl(7); rd = pl(8); Vf = pl(9);
[ic, vo] = cap(iL, vC, ip, R, rC);
rho = rL + d*Ron + (1 - d)*rd;
diL = (d*Vin - rho*iL - (1 - d)*Vf - vo)/L;
if iL <= 0 && diL < 0, diL = 0; end
dvC = ic/C;
pin = d*Vin*iL;
pu = vo*vo/R + vo*ip;
pls = rho*iL*iL + (1 - d)*Vf*iL + rC*ic*ic;
end

function out = sim_avg(P, p, tr, lv, ld, iL0, vC0, T, nsub)
% Modèle moyen, Eq. (6), garde i_L >= 0 (DCM moyen), RK4, commande à Ts.
% ld = [t_on t_off amplitude]. out : t, r, vo, iL, d, dpre, diag
pl = plant_vec(P); R = pl(4); rC = pl(6);
Ts = p(33); n = round(T/Ts); h = Ts/nsub;
iL = iL0; vC = vC0;
[~, vo] = cap(iL, vC, 0, R, rC);
s = ctrl_init(p, ref(0, tr, lv), vo, iL, vC);
tv = zeros(n+1,1); rv = tv; vv = tv; iv = tv; dv = tv; dp = tv;
Ein = 0; Eu = 0; El = 0; nclip = 0; nsat = 0; ipk = iL;
out.ok = true;
for k = 0:n-1
    t = k*Ts; r = ref(t, tr, lv);
    ip = 0; if t >= ld(1) && t < ld(2), ip = ld(3); end
    [ic, vo] = cap(iL, vC, ip, R, rC);
    [s, dnow, dsat, dt_] = ctrl_step(p, s, r, vo, iL, ic);
    tv(k+1) = t; rv(k+1) = r; vv(k+1) = vo; iv(k+1) = iL; dv(k+1) = dnow; dp(k+1) = dt_;
    if dsat ~= dt_, nsat = nsat + 1; end
    for j = 1:nsub
        [a1, b1, p1, u1, l1] = avg_f(iL, vC, dnow, ip, pl);
        [a2, b2, p2, u2, l2] = avg_f(iL + 0.5*h*a1, vC + 0.5*h*b1, dnow, ip, pl);
        [a3, b3, p3, u3, l3] = avg_f(iL + 0.5*h*a2, vC + 0.5*h*b2, dnow, ip, pl);
        [a4, b4, p4, u4, l4] = avg_f(iL + h*a3, vC + h*b3, dnow, ip, pl);
        iL = iL + h/6*(a1 + 2*a2 + 2*a3 + a4);
        vC = vC + h/6*(b1 + 2*b2 + 2*b3 + b4);
        Ein = Ein + h/6*(p1 + 2*p2 + 2*p3 + p4);
        Eu = Eu + h/6*(u1 + 2*u2 + 2*u3 + u4);
        El = El + h/6*(l1 + 2*l2 + 2*l3 + l4);
        if iL < 0, iL = 0; nclip = nclip + 1; end
        if iL > ipk, ipk = iL; end
    end
    if ~(abs(vC) < 1e4 && abs(iL) < 1e4)
        out.ok = false; break;
    end
end
if out.ok
    t = n*Ts; r = ref(t, tr, lv);
    ip = 0; if t >= ld(1) && t < ld(2), ip = ld(3); end
    [~, vo] = cap(iL, vC, ip, R, rC);
    tv(n+1) = t; rv(n+1) = r; vv(n+1) = vo; iv(n+1) = iL; dv(n+1) = s(5); dp(n+1) = s(6);
end
out.t = tv; out.r = rv; out.vo = vv; out.iL = iv; out.d = dv; out.dpre = dp;
out.Ein = Ein; out.Eu = Eu; out.El = El; out.clip = nclip/(n*nsub);
out.sat = nsat/n; out.ipk = ipk;
end

function [diL, dvC, pin, pu, pls, mode] = sw_f(iL, vC, q, ip, pl, fm)
% Réplique de l'étage de puissance du .slx (Fig. 2) ; mode 0 on, 1 diode,
% 2 DCM (branche de fuite quasi statique), 3 diode interne du MOSFET
Vin = pl(1); L = pl(2); C = pl(3); R = pl(4); rL = pl(5); rC = pl(6);
Ron = pl(7); rd = pl(8); Vf = pl(9); Rs = pl(10); Rb = pl(11);
[ic, vo] = cap(iL, vC, ip, R, rC);
ith = (Vin + Vf)/Rs;
if fm >= 0
    mode = fm;
elseif q == 1
    mode = 0;
elseif iL < 0
    mode = 3;
elseif iL > ith
    mode = 1;
else
    mode = 2;
end
vd = 0; idd = 0;
switch mode
    case 0
        gs = 1/Ron + 1/Rs;
        if iL < 0, gs = gs + 1/Rb; end
        vx = Vin - iL/gs; isw = iL;
    case 3
        vx = Vin - iL/(1/Rb + 1/Rs); isw = iL;
    case 1
        vx = (Vin/Rs - Vf/rd - iL)/(1/Rs + 1/rd);
        idd = (-vx - Vf)/rd; vd = -vx;
        isw = (Vin - vx)/Rs;
    otherwise
        vx = Vin - Rs*iL; isw = iL;
end
diL = (vx - rL*iL - vo)/L;
if mode == 2, diL = 0; end
dvC = ic/C;
pin = Vin*isw;
pu = vo*vo/R + vo*ip;
pls = rL*iL*iL + rC*ic*ic + (Vin - vx)*isw + vd*idd;
end

function out = sim_sw(P, p, tr, lv, ld, iL0, vC0, T, h, win)
% Réplique commutée : fronts MLI exacts, diode bloquante (DCM naturel),
% mode figé dans chaque pas RK4, bilan d'énergie fermé.
pl = plant_vec(P);
Vin = pl(1); L = pl(2); C = pl(3); R = pl(4); rL = pl(5); rC = pl(6); Vf = pl(9); Rs = pl(10);
Ts = p(33); Tsw = 1/P.fs; n = round(T/Ts);
if nargin < 10, win = zeros(0, 2); end
nw = size(win, 1);
iL = iL0; vC = vC0;
[~, vo] = cap(iL, vC, 0, R, rC);
s = ctrl_init(p, ref(0, tr, lv), vo, iL, vC);
tv = zeros(n+1,1); rv = tv; vv = tv; iv = tv; dv = tv; dp = tv; fdcm = tv;
Ein = 0; Eu = 0; El = 0;
Es0 = 0.5*L*iL^2 + 0.5*C*vC^2;
wacc = zeros(nw, 5); wacc(:, 4) = NaN;
ipk = iL; vmax = vo; t_dcm = 0; t_rev = 0; nsat = 0; iae = 0;
ep = 1e-13;
for k = 0:n-1
    t = k*Ts; r = ref(t, tr, lv);
    ip = 0; if t >= ld(1) && t < ld(2), ip = ld(3); end
    [ic, vo] = cap(iL, vC, ip, R, rC);
    [s, dnow, dsat, dt_] = ctrl_step(p, s, r, vo, iL, ic);
    tv(k+1) = t; rv(k+1) = r; vv(k+1) = vo; iv(k+1) = iL; dv(k+1) = dnow; dp(k+1) = dt_;
    if dsat ~= dt_, nsat = nsat + 1; end
    for w = 1:nw
        if abs(t - win(w,1)) < 0.5*Ts
            wacc(w, 4) = 0.5*L*iL^2 + 0.5*C*vC^2;
            wacc(w, 1:3) = [Ein Eu El];
        end
    end
    ton = dnow*Tsw; tdk = t_dcm; tl = t; tstop = t + Ts;
    while tl < tstop - ep
        nper = floor((tl + ep)/Tsw);
        ph = tl - nper*Tsw;
        if ph < ton - ep
            q = 1; tedge = nper*Tsw + ton;
        else
            q = 0; tedge = (nper + 1)*Tsw;
        end
        if tedge <= tl + ep, tedge = tl + Tsw; end
        dt = min([tstop - tl, tedge - tl, h]);
        if dt <= ep, dt = tstop - tl; end
        e0 = abs(r - vo);
        [a1, b1, p1, u1, l1, m1] = sw_f(iL, vC, q, ip, pl, -1);
        if m1 == 2
            [~, b2, p2, u2, l2] = sw_f(iL, vC + 0.5*dt*b1, q, ip, pl, m1);
            [~, b3, p3, u3, l3] = sw_f(iL, vC + 0.5*dt*b2, q, ip, pl, m1);
            [~, b4, p4, u4, l4] = sw_f(iL, vC + dt*b3, q, ip, pl, m1);
            vC = vC + dt/6*(b1 + 2*b2 + 2*b3 + b4);
            Ein = Ein + dt/6*(p1 + 2*p2 + 2*p3 + p4);
            Eu = Eu + dt/6*(u1 + 2*u2 + 2*u3 + u4);
            El = El + dt/6*(l1 + 2*l2 + 2*l3 + l4);
            [~, vo] = cap(iL, vC, ip, R, rC);
            inew = (Vin - vo)/(Rs + rL);
            El = El + 0.5*L*(iL^2 - inew^2);
            iL = inew; t_dcm = t_dcm + dt;
        else
            [a2, b2, p2, u2, l2] = sw_f(iL + 0.5*dt*a1, vC + 0.5*dt*b1, q, ip, pl, m1);
            [a3, b3, p3, u3, l3] = sw_f(iL + 0.5*dt*a2, vC + 0.5*dt*b2, q, ip, pl, m1);
            [a4, b4, p4, u4, l4] = sw_f(iL + dt*a3, vC + dt*b3, q, ip, pl, m1);
            iL = iL + dt/6*(a1 + 2*a2 + 2*a3 + a4);
            vC = vC + dt/6*(b1 + 2*b2 + 2*b3 + b4);
            Ein = Ein + dt/6*(p1 + 2*p2 + 2*p3 + p4);
            Eu = Eu + dt/6*(u1 + 2*u2 + 2*u3 + u4);
            El = El + dt/6*(l1 + 2*l2 + 2*l3 + l4);
            if q == 0 && ((m1 == 1 && iL < (Vin + Vf)/Rs) || (m1 == 3 && iL > 0))
                [~, vo] = cap(iL, vC, ip, R, rC);
                inew = (Vin - vo)/(Rs + rL);
                El = El + 0.5*L*(iL^2 - inew^2);
                iL = inew;
            end
            if iL < 0, t_rev = t_rev + dt; end
        end
        [~, vo] = cap(iL, vC, ip, R, rC);
        iae = iae + 0.5*(e0 + abs(r - vo))*dt;
        if iL > ipk, ipk = iL; end
        if vo > vmax, vmax = vo; end
        tl = tl + dt;
    end
    fdcm(k+1) = (t_dcm - tdk)/Ts;
    for w = 1:nw
        if abs(tstop - win(w,2)) < 0.5*Ts
            Es = 0.5*L*iL^2 + 0.5*C*vC^2;
            wacc(w, 1:3) = [Ein Eu El] - wacc(w, 1:3);
            wacc(w, 5) = Es - wacc(w, 4);
        end
    end
    if ~(abs(vC) < 1e4 && abs(iL) < 1e4), break; end
end
t = n*Ts; r = ref(t, tr, lv);
[~, vo] = cap(iL, vC, 0, R, rC);
tv(n+1) = t; rv(n+1) = r; vv(n+1) = vo; iv(n+1) = iL; dv(n+1) = s(5); dp(n+1) = s(6);
dEs = 0.5*L*iL^2 + 0.5*C*vC^2 - Es0;
out.t = tv; out.r = rv; out.vo = vv; out.iL = iv; out.d = dv; out.dpre = dp; out.fdcm = fdcm;
out.Ein = Ein; out.Eu = Eu; out.Eloss = El; out.dEs = dEs; out.res = Ein - Eu - El - dEs;
out.res_pct = 100*out.res/Ein; out.ipk = ipk; out.vmax = vmax; out.dcm = t_dcm/T;
out.sat = nsat/n; out.trev = t_rev; out.IAE = iae; out.eta = Eu/Ein;
out.res_w = zeros(1, nw); out.eta_w = zeros(1, nw);
for w = 1:nw
    out.res_w(w) = 100*(wacc(w,1) - wacc(w,2) - wacc(w,3) - wacc(w,5))/wacc(w,1);
    out.eta_w(w) = wacc(w,2)/wacc(w,1);
end
end

function out = mission_sw(P, p)
% mission de validation complète sur la réplique commutée (pas 0.2 us)
m = P.mission;
out = sim_sw(P, p, m.t_steps, m.levels, [m.load m.load_amp], m.levels(1)/P.R, m.levels(1), ...
             m.T, m.h, m.win);
end

%% ===================================================== F. évaluateur
function m = scenario_metrics(P, p, sc)
IL0 = equilibrium(P, sc.lv(1));
o = sim_avg(P, p, sc.tr, sc.lv, sc.load, IL0, sc.lv(1), sc.T, P.search.nsub);
if ~o.ok, m = []; return; end
tv = o.t; e = o.r - o.vo; ae = abs(e); T = sc.T;
m.IAE = trapz(tv, ae); m.ISE = trapz(tv, e.^2); m.ITAEg = trapz(tv, tv.*ae);
rs = [sc.resets T];
m.ITAEev = 0; ess = 0;
for k = 1:numel(rs) - 1
    a = rs(k); b = rs(k+1);
    sel = tv >= a & tv <= b;
    m.ITAEev = m.ITAEev + trapz(tv(sel), (tv(sel) - a).*ae(sel));
    w = tv >= b - P.spec.ESS_WIN & tv < b;
    ess = max(ess, abs(mean(e(w))));
end
os = -Inf;
for k = 1:numel(sc.tr)
    sg = sign(sc.lv(k+1) - sc.lv(k));
    if k < numel(sc.tr), nxt = sc.tr(k+1); else, nxt = T; end
    sel = tv >= sc.tr(k) & tv < nxt;
    os = max(os, max(sg*(o.vo(sel) - sc.lv(k+1))));
end
ld = -Inf;
if sc.load(3) > 0
    cand = rs(rs > sc.load(2) + 1e-12);
    nxt = min([cand T]);
    sel = tv >= sc.load(1) & tv < nxt;
    ld = max(ae(sel));
end
m.RMSE = sqrt(m.ISE/T); m.MAE = m.IAE/T; m.ILpk = o.ipk; m.ess = ess;
m.Eu = o.Eu; m.El = o.El; m.LI = o.El/o.Eu; m.clip = o.clip; m.sat = o.sat;
m.OS = os; m.LD = ld;
end

function [dpp, ess, ipk] = sw_screen(P, p)
% écran commuté : Vmid -> Vlow -> Vhigh -> Vmid en 34 ms (RK4 1 us)
S = P.screen;
IL = equilibrium(P, S.lv(1));
o = sim_sw(P, p, S.tr, S.lv, [0 0 0], IL, S.lv(1), S.T, P.search.h_screen, zeros(0, 2));
if ~all(isfinite(o.vo)), dpp = Inf; ess = Inf; ipk = Inf; return; end
dpp = 0; ess = 0;
for w = 1:size(S.win, 1)
    sel = o.t >= S.win(w,1) & o.t < S.win(w,2);
    dpp = max(dpp, max(o.d(sel)) - min(o.d(sel)));
    ess = max(ess, abs(mean(o.r(sel) - o.vo(sel))));
end
ipk = o.ipk;
end

function E = evaluator(P, dec, p_ref)
% E.fn(x) renvoie un enregistrement de Deb ; normalisation par p_ref (ancre xi0)
if nargin < 3, p_ref = dec_cmo(P, P.lqi.xi0); end
E.P = P; E.dec = dec;
E.ref = cell(1, numel(P.scen));
for k = 1:numel(P.scen)
    E.ref{k} = scenario_metrics(P, p_ref, P.scen(k));
end
E.fn = @(x) score(E, x);
end

function rec = score(E, x)
P = E.P; S = P.spec;
p = E.dec(x);
crit = {'LI','IAE','ISE','ITAEg','ITAEev','RMSE','MAE'};
ns = numel(P.scen);
rec = struct('x', x(:)', 'u', [], 'J', 1e3, 'g', 1e3*ones(1, 33), 'viol', 33e3, ...
             'feasible', false, 'dpp', Inf, 'ess_sw', Inf, 'ipk_sw', Inf, ...
             'ratios', [], 'n_eval', 0, 'stage', '');
ms = cell(1, ns);
for k = 1:ns
    ms{k} = scenario_metrics(P, p, P.scen(k));
    if isempty(ms{k}), return; end
end
r = zeros(ns, 7); g = zeros(1, 33); j = 0;
for k = 1:ns
    m = ms{k}; m0 = E.ref{k};
    for c = 1:7, r(k, c) = m.(crit{c})/m0.(crit{c}); end
    if isfinite(m.OS), g(j+1) = m.OS/S.OS_LIM - 1; else, g(j+1) = -1; end
    if isfinite(m.LD), g(j+2) = m.LD/S.LOAD_LIM - 1; else, g(j+2) = -1; end
    g(j+3) = (m.ILpk + S.HALF_RIPPLE)/S.I_LIM - 1;
    g(j+4) = m.ess/S.ESS_TOL - 1;
    g(j+5) = (m.LI/m0.LI - 1)/S.LI_TOL - 1;
    j = j + 5;
end
[dpp, ess_sw, ipk_sw] = sw_screen(P, p);
g(31) = dpp/S.DPP_LIM - 1; g(32) = ess_sw/S.ESS_TOL - 1; g(33) = ipk_sw/S.I_LIM - 1;
rec.J = mean(r(:)); rec.g = g; rec.viol = sum(max(g, 0));
rec.feasible = rec.viol <= S.FEAS_TOL;
rec.dpp = dpp; rec.ess_sw = ess_sw; rec.ipk_sw = ipk_sw; rec.ratios = r;
end

%% ================================================ G. algorithmes 1-7
function tf = better(a, b)
if isempty(b), tf = true; return; end
if a.feasible ~= b.feasible, tf = a.feasible; return; end
if a.feasible, tf = a.J < b.J; else, tf = a.viol < b.viol; end
end

function b = best_of(recs)
b = [];
for k = 1:numel(recs)
    if better(recs{k}, b), b = recs{k}; end
end
end

function o = order(recs)
n = numel(recs); key = zeros(n, 2);
for k = 1:n
    if recs{k}.feasible, key(k, :) = [0 recs{k}.J]; else, key(k, :) = [1 recs{k}.viol]; end
end
[~, o] = sortrows(key);
o = o(:)';
end

function C = counter(fn, lb, ub, budget, stage)
% compteur d'appels partagé (containers.Map est un objet « handle »)
C = containers.Map();
C('n') = 0; C('hist') = {};
C('meta') = struct('fn', fn, 'lb', lb(:)', 'ub', ub(:)', 'budget', budget, 'stage', stage);
end

function rec = cev(C, u)
M = C('meta');
if C('n') >= M.budget, error('cmo:budget', 'budget épuisé'); end
u = min(max(u(:)', 0), 1);
rec = M.fn(M.lb + u.*(M.ub - M.lb));
rec.u = u;
n = C('n') + 1; C('n') = n;
rec.n_eval = n; rec.stage = M.stage;
h = C('hist'); h{end+1} = rec; C('hist') = h;
end

function t = ctau(C)
M = C('meta'); t = min(1, C('n')/M.budget);
end

function seed_rng(seed)
if exist('OCTAVE_VERSION', 'builtin')
    rand('twister', seed); randn('twister', seed); %#ok<RAND>
else
    rng(seed, 'twister');
end
end

function st = rng_state(st)
% sans argument : lit l'état ; avec argument : le restaure
if exist('OCTAVE_VERSION', 'builtin')
    if nargin, rand('twister', st{1}); randn('twister', st{2}); else, st = {rand('twister'), randn('twister')}; end %#ok<RAND>
else
    if nargin, rng(st); else, st = rng; end
end
end

function u = reflect(u)
u(u < 0) = -u(u < 0);
u(u > 1) = 2 - u(u > 1);
u = min(max(u, 0), 1);
end

function X = init_pop(u0, d, NPOP)
X = rand(NPOP, d);
nloc = round(0.6*NPOP);
for i = 2:nloc
    z = randn(1, d);
    X(i, :) = u0 + 0.05*rand*z/norm(z);
end
X(1, :) = u0;
X = reflect(X);
end

function pop = eval_pop(C, X)
pop = cell(1, size(X, 1));
for i = 1:size(X, 1), pop{i} = cev(C, X(i, :)); end
end

function alg_pso(C, u0, NPOP)
d = numel(u0);
X = init_pop(u0, d, NPOP); V = zeros(size(X));
Pb = eval_pop(C, X); Px = X; G = best_of(Pb);
while true
    w = 0.9 - 0.5*ctau(C);
    for i = 1:NPOP
        V(i,:) = w*V(i,:) + 1.7*rand(1,d).*(Px(i,:) - X(i,:)) + 1.7*rand(1,d).*(G.u - X(i,:));
        V(i,:) = min(max(V(i,:), -0.2), 0.2);
        X(i,:) = reflect(X(i,:) + V(i,:));
        r = cev(C, X(i,:));
        if better(r, Pb{i}), Pb{i} = r; Px(i,:) = X(i,:); end
        if better(r, G), G = r; end
    end
end
end

function alg_ga(C, u0, NPOP)
d = numel(u0);
pop = eval_pop(C, init_pop(u0, d, NPOP));
while true
    kids = {};
    for q = 1:NPOP/2
        a = ga_tour(pop, NPOP); b = ga_tour(pop, NPOP);
        a = a.u; b = b.u;
        uu = rand(1, d);
        beta = (2*uu).^(1/16);
        beta(uu > 0.5) = (1./(2*(1 - uu(uu > 0.5)))).^(1/16);
        ch = {0.5*((1 + beta).*a + (1 - beta).*b), 0.5*((1 - beta).*a + (1 + beta).*b)};
        for c = 1:2
            y = ch{c};
            for j = 1:d
                if rand < 1/d
                    qq = rand;
                    if qq < 0.5, dl = (2*qq)^(1/21) - 1; else, dl = 1 - (2*(1 - qq))^(1/21); end
                    y(j) = y(j) + dl;
                end
            end
            kids{end+1} = cev(C, reflect(y)); %#ok<AGROW>
        end
    end
    merged = [pop kids];
    o = order(merged);
    pop = merged(o(1:NPOP));
end
end

function r = ga_tour(pop, NPOP)
i = randi(NPOP); j = randi(NPOP);
if better(pop{i}, pop{j}), r = pop{i}; else, r = pop{j}; end
end

function alg_abc(C, u0, NPOP, adaptive)
d = numel(u0);
pop = eval_pop(C, init_pop(u0, d, NPOP));
trials = zeros(1, NPOP);
limit = max(5, round(0.5*NPOP*d));
while true
    tau = ctau(C); sc = 1;
    if adaptive
        U = cell2mat(cellfun(@(r) r.u, pop, 'UniformOutput', false)');
        Dk = mean(std(U, 1, 1));
        sc = min(max(0.5 + 0.15/max(Dk, 1e-3), 0.5), 2);
    end
    for i = 1:NPOP
        y = abc_move(pop, i, sc, NPOP, d);
        if adaptive && rand < 0.25*(1 - tau)
            j2 = randi(d); y(j2) = y(j2) + (rand*0.2 - 0.1);
        end
        r = cev(C, reflect(y));
        if better(r, pop{i}), pop{i} = r; trials(i) = 0; else, trials(i) = trials(i) + 1; end
    end
    o = order(pop);
    rk = zeros(1, NPOP); rk(o) = 0:NPOP-1;
    pr = (NPOP - rk)/sum(NPOP - rk);
    cp = cumsum(pr);
    for q = 1:NPOP
        i = find(rand <= cp, 1); if isempty(i), i = NPOP; end
        r = cev(C, reflect(abc_move(pop, i, sc, NPOP, d)));
        if better(r, pop{i}), pop{i} = r; trials(i) = 0; else, trials(i) = trials(i) + 1; end
    end
    [tm, i] = max(trials);
    if tm >= limit, pop{i} = cev(C, rand(1, d)); trials(i) = 0; end
end
end

function y = abc_move(pop, i, sc, NPOP, d)
k = randi(NPOP - 1); if k >= i, k = k + 1; end
j = randi(d);
y = pop{i}.u;
y(j) = y(j) + sc*(2*rand - 1)*(y(j) - pop{k}.u(j));
end

function L = leaders(pop)
o = order(pop);
L = {pop{o(1)}.u, pop{o(2)}.u, pop{o(3)}.u};
end

function acc = gwo_move(x, L, a, d)
acc = zeros(1, d);
for m = 1:3
    A = 2*a*rand(1, d) - a; Cm = 2*rand(1, d);
    acc = acc + L{m} - A.*abs(Cm.*L{m} - x);
end
acc = acc/3;
end

function alg_igwo(C, u0, NPOP)
d = numel(u0);
pop = eval_pop(C, init_pop(u0, d, NPOP));
while true
    tau = ctau(C); a = 2*(1 - tau^2);
    L = leaders(pop); cand = cell(1, NPOP);
    for i = 1:NPOP
        y = gwo_move(pop{i}.u, L, a, d);
        if rand < 0.15*(1 - tau) + 0.02
            j = randi(d); y(j) = y(j) + 0.1*randn;
        end
        cand{i} = cev(C, reflect(y));
    end
    merged = [pop cand];
    o = order(merged);
    pop = merged(o(1:NPOP));
end
end

function alg_igwo_dlh(C, u0, NPOP)
d = numel(u0); M = C('meta');
pop = eval_pop(C, init_pop(u0, d, NPOP));
K = ceil((M.budget - NPOP)/(2*NPOP)); k = 0;
while true
    a = max(0, 2 - 2*k/K);
    L = leaders(pop);
    for i = 1:NPOP
        r1 = cev(C, min(max(gwo_move(pop{i}.u, L, a, d), 0), 1));
        Ri = 0.5*norm(pop{i}.u - L{1});
        z = randn(1, d);
        r2 = cev(C, min(max(pop{i}.u + Ri*rand*z/norm(z), 0), 1));
        if better(r1, r2), pop{i} = r1; else, pop{i} = r2; end
    end
    k = k + 1;
end
end

function K = matern52(A, B, ell)
D2 = max(sum(A.^2, 2) + sum(B.^2, 2)' - 2*A*B', 0);
r = sqrt(D2)/ell; s5 = sqrt(5)*r;
K = (1 + s5 + 5*r.^2/3).*exp(-s5);
end

function gp = gp_fit(U, Y)
% GP à noyau Matérn 5/2 partagé par J et les 33 contraintes
mu = mean(Y, 1); sd = std(Y, 1, 1); sd(sd < 1e-12) = 1;
Z = (Y - mu)./sd; n = size(U, 1);
best = -Inf;
for ell = [0.05 0.1 0.2 0.4 0.8]
    [Lc, fl] = chol(matern52(U, U, ell) + 1e-6*eye(n), 'lower');
    if fl, continue; end
    a = Lc'\(Lc\Z(:,1));
    ll = -0.5*Z(:,1)'*a - sum(log(diag(Lc)));
    if ll > best, best = ll; gp.ell = ell; gp.Lc = Lc; end
end
gp.U = U; gp.mu = mu; gp.sd = sd;
gp.alpha = gp.Lc'\(gp.Lc\Z);
end

function [m, s] = gp_pred(gp, V)
Ks = matern52(V, gp.U, gp.ell);
m = (Ks*gp.alpha).*gp.sd + gp.mu;
v = gp.Lc\Ks';
var = max(1 - sum(v.^2, 1)', 1e-12);
s = sqrt(var)*gp.sd;
end

function alg_refine_bo(C, hist)
% Algorithme 7 : GP + EI pondérée par la faisabilité + région de confiance
M = C('meta'); d = numel(hist{1}.u);
allr = hist;
Phi = @(z) 0.5*erfc(-z/sqrt(2)); phi = @(z) exp(-z.^2/2)/sqrt(2*pi);
for k = 0:M.budget-1
    n = numel(allr);
    U = zeros(n, d); Y = zeros(n, 34);
    for i = 1:n
        U(i,:) = allr{i}.u;
        Y(i,:) = [min(allr{i}.J, 5), min(max(allr{i}.g, -5), 5)];
    end
    gp = gp_fit(U, Y);
    rad = 0.04*(0.003/0.04)^(k/max(M.budget - 1, 1));
    o = order(allr); o = o(1:min(10, n));
    Cn = U(o, :);
    z = randn(1000, d);
    z = z.*((rand(1000, 1).^(1/d))./sqrt(sum(z.^2, 2)));
    pool = min(max(Cn(randi(numel(o), 1000, 1), :) + rad*z, 0), 1);
    [m, s] = gp_pred(gp, pool);
    mJ = m(:,1); sJ = s(:,1); mg = m(:,2:end); sg = s(:,2:end);
    Pf = prod(Phi(-mg./sg), 2);
    inc = best_of(allr);
    if inc.feasible
        zz = (inc.J - mJ)./sJ;
        EI = (inc.J - mJ).*Phi(zz) + sJ.*phi(zz);
        [~, i] = max(Pf.*(EI + 0.02*sJ) + 1e-4*Pf);
    else                               % mode de restauration prédéclaré
        [~, i] = min(sum(max(mg, 0), 2) - 0.1*mean(sg, 2));
    end
    allr{end+1} = cev(C, pool(i, :)); %#ok<AGROW>
end
end

function hist = run_method(P, method, fn, lb, ub, x0, seed)
% Une exécution à graine appariée : B_global appels globaux + B_local refineBO
lb = lb(:)'; ub = ub(:)';
u0 = (x0(:)' - lb)./(ub - lb);
seed_rng(seed);
NPOP = P.search.npop;
Cg = counter(fn, lb, ub, P.search.B_global, 'global');
try
    switch method
        case 'PSO',       alg_pso(Cg, u0, NPOP);
        case 'GA',        alg_ga(Cg, u0, NPOP);
        case 'ABC',       alg_abc(Cg, u0, NPOP, false);
        case 'pAEABC',    alg_abc(Cg, u0, NPOP, true);
        case 'pIGWO',     alg_igwo(Cg, u0, NPOP);
        case 'pIGWO-DLH', alg_igwo_dlh(Cg, u0, NPOP);
        otherwise, error('cmo:method', 'méthode inconnue : %s', method);
    end
catch err
    if ~strcmp(err.identifier, 'cmo:budget'), rethrow(err); end
end
hg = Cg('hist');
Cl = counter(fn, lb, ub, P.search.B_local, 'local');
if P.search.B_local > 0
    try
        alg_refine_bo(Cl, hg);
    catch err
        if ~strcmp(err.identifier, 'cmo:budget'), rethrow(err); end
    end
end
hl = Cl('hist');
for i = 1:numel(hl), hl{i}.n_eval = hl{i}.n_eval + P.search.B_global; end
hist = [hg hl];
end

%% ======================================= H. rejeu sur le modèle Simulink
function [m, so] = sim_simulink(P, p, folder)
% Rejoue le correcteur p sur le .slx de référence avec les paramètres de P.
% Le .slx n'est jamais modifié : tout passe par Simulink.SimulationInput.
if nargin < 3, folder = fileparts(mfilename('fullpath')); end
mdl = P.slx;
load_system(fullfile(folder, [mdl '.slx']));
[blk, vars] = model_settings(P, p);
in = Simulink.SimulationInput(mdl);
for k = 1:size(blk, 1)
    in = in.setBlockParameter([mdl '/' blk{k,1}], blk{k,2}, blk{k,3});
end
for k = 1:size(vars, 1)
    in = in.setVariable(vars{k,1}, vars{k,2}, 'Workspace', mdl);
end
in = in.setModelParameter('StopTime', num2str(P.mission.T, 17));
so = sim(in);
m = mission_metrics(P, so);
end

function configure_model(P, p, folder)
% Ouvre le .slx et y applique P et p EN MÉMOIRE (fichier non sauvegardé) :
% on peut alors cliquer sur Run et regarder les Scopes.
if nargin < 3, folder = fileparts(mfilename('fullpath')); end
mdl = P.slx;
open_system(fullfile(folder, [mdl '.slx']));
[blk, vars] = model_settings(P, p);
for k = 1:size(blk, 1)
    set_param([mdl '/' blk{k,1}], blk{k,2}, blk{k,3});
end
mws = get_param(mdl, 'ModelWorkspace');
for k = 1:size(vars, 1)
    assignin(mws, vars{k,1}, vars{k,2});
end
set_param(mdl, 'StopTime', num2str(P.mission.T, 17));
fprintf('%s configuré (non sauvegardé) : cliquer sur Run.\n', mdl);
end

function [blk, vars] = model_settings(P, p)
% paramètres de blocs et variables de l'espace de travail du modèle
ms = P.mission; V0 = ms.levels(1);
n2s = @(v) num2str(v, 17);
blk = {P.ctrlblk, 'FunctionName', 'cmo_sfun';
       P.ctrlblk, 'Parameters', mat2str(p, 17);
       'L', 'Inductance', n2s(P.L);
       'L', 'InitialCurrent', n2s(V0/P.R);
       'C', 'Capacitance', n2s(P.C);
       'C', 'InitialVoltage', n2s(V0);
       'R', 'Resistance', n2s(P.R);
       'rl', 'Resistance', n2s(P.rL);
       'R1', 'Resistance', n2s(P.rC);
       'Vin constant 48 V', 'Before', n2s(P.Vin);
       'Vin constant 48 V', 'After', n2s(P.Vin);
       'IGBT/Diode', 'Ron', n2s(P.Ron);
       'IGBT/Diode', 'Rd', n2s(P.Rbody);
       'IGBT/Diode', 'Rs', n2s(P.Rs);
       'Diode1', 'Ron', n2s(P.rd);
       'Diode1', 'Vf', n2s(P.Vf);
       'Resistor current v_o over R', 'Gain', n2s(1/P.R);
       P.pwmblk, 'Fsw', n2s(P.fs)};
tr = [0 ms.t_steps ms.T]; lvv = [ms.levels ms.levels(end)];
tl = [0 ms.load ms.T]; il = [0 ms.load_amp 0 0];
vars = {'Vin', P.Vin; 'L', P.L; 'C', P.C; 'R', P.R; 'rL', P.rL; 'rC', P.rC; ...
        'Ron', P.Ron; 'Rs', P.Rs; 'rd', P.rd; 'Vf', P.Vf; 'fs', P.fs; ...
        'd_min', P.dmin; 'd_max', P.dmax; 'Tsim', ms.T; 'CTRL_PARAM_VECTOR', p; ...
        'Vref_ts', timeseries(lvv(:), tr(:)); 'Iload_ts', timeseries(il(:), tl(:))};
end

function m = mission_metrics(P, so)
% indices de la mission (mêmes définitions que la réplique et que l'article)
ms = P.mission;
t = so.V_o_log.Time; vo = so.V_o_log.Data(:);
iL = rs1(so.i_L_log, t); iC = rs1(so.i_C_log, t); iin = rs1(so.i_in_log, t);
ip = rs1(so.i_p_log, t);
iSW = rs2(so.igbt_m_log, t, 1); vSW = rs2(so.igbt_m_log, t, 2);
iD = rs2(so.diode_m_log, t, 1); vD = rs2(so.diode_m_log, t, 2);
r = ms.levels(1)*(t < ms.t_steps(1)) + ms.levels(2)*(t >= ms.t_steps(1) & t < ms.t_steps(2)) ...
    + ms.levels(3)*(t >= ms.t_steps(2));
vC = vo - P.rC*iC;
Pin = P.Vin*iin; Pu = vo.^2/P.R + vo.*ip;
Pl = P.rL*iL.^2 + P.rC*iC.^2 + vSW.*iSW + vD.*iD + vSW.^2/P.Rs;
Es = 0.5*P.L*iL.^2 + 0.5*P.C*vC.^2;
Ein = trapz(t, Pin); Eu = trapz(t, Pu); El = trapz(t, Pl);
m.IAE = trapz(t, abs(r - vo));
m.ipk = max(iL); m.vmax = max(vo);
m.os_high = max(vo(t >= ms.t_steps(2))) - ms.levels(3);
for w = 1:3
    sel = t >= ms.win(w,1) & t <= ms.win(w,2);
    ts = t(sel);
    m.(sprintf('ess%d', w)) = trapz(ts, r(sel) - vo(sel))/(ts(end) - ts(1));
    m.(sprintf('eta%d', w)) = trapz(ts, Pu(sel))/trapz(ts, Pin(sel));
end
dsat = so.d_sat_log.Data(:); dpre = so.d_pre_log.Data(:);
m.sat = mean(abs(dsat - dpre) > 1e-12);
m.eta = Eu/Ein; m.LI = El/Eu; m.Eu = Eu; m.Eloss = El;
m.res_pct = 100*(Ein - Eu - El - (Es(end) - Es(1)))/Ein;
end

function y = rs1(ts, t)
[tt, i] = unique(ts.Time, 'last'); d = ts.Data(:); d = d(i);
y = interp1(tt, d, t, 'linear', 'extrap');
end

function y = rs2(ts, t, col)
[tt, i] = unique(ts.Time, 'last'); d = squeeze(ts.Data);
if size(d, 1) ~= numel(ts.Time), d = d.'; end
y = interp1(tt, d(i, col), t, 'linear', 'extrap');
end
