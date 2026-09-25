function fopidf_sfun(block)
% FOPIDF_SFUN  Correcteur FO-PIDF de la campagne v6.1, en S-function niveau 2.
%
% Remplace buck_comparator_sfun dans le bloc "Common sampled controller" du
% modele Simulink, avec la MEME interface :
%   entrees : 1 Vref, 2 V_o, 3 I_L, 4 I_c, 5 V_in   (seuls Vref et V_o servent)
%   sorties : 1 d_sat (vers la MLI), 2 d_pre (avant saturation), 3 lambda
%
% Traduction ligne a ligne de isafo/_kernel.py : Oustaloup demi-ordre 3
% (7 cellules) sur [1 ; 1e5] rad/s discretise exactement par BOZ, integrale
% s^-lambda = (1/s) * s^(1-lambda), derivee s^mu filtree par wf/(s+wf),
% anti-emballement par recalcul, feedforward nominal, RETARD D'UN ECHANTILLON.
%
% Parametre unique (vecteur, 19 elements) :
%   [kp ki kd lambda mu wf b c Kb Ts Vn Dn dmin dmax R rL rd Ron Vin]
% R, rL, rd, Ron, Vin sont les valeurs NOMINALES du feedforward, figees meme
% quand le plant varie (par.1 du cadrage).

setup(block);
end

function setup(block)
block.NumInputPorts  = 5;
block.NumOutputPorts = 3;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
for k = 1:5
    block.InputPort(k).Dimensions = 1;
    block.InputPort(k).DirectFeedthrough = false;  % retard d'un echantillon
end
for k = 1:3
    block.OutputPort(k).Dimensions = 1;
end
block.NumDialogPrms = 1;
P = block.DialogPrm(1).Data;
block.SampleTimes = [P(10) 0];
block.SimStateCompliance = 'DefaultSimState';
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitConditions);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Update',  @Update);
end

function DoPostPropSetup(block)
block.NumDworks = 5;
names = {'wi', 'wd', 'xs', 'dout', 'coef'};
dims  = [7, 7, 2, 2, 4*7 + 4*7 + 2];
for k = 1:5
    block.Dwork(k).Name = names{k};
    block.Dwork(k).Dimensions = dims(k);
    block.Dwork(k).DatatypeID = 0;
    block.Dwork(k).Complexity = 'Real';
    block.Dwork(k).UsedAsDiscState = true;
end
end

function [z, p, g] = oustaloup(gam)
% Meme formule que isafo/oustaloup.py.
N = 3; wb = 1; wh = 1e5; n = 2*N + 1; k = (-N:N)';
gam = min(max(gam, 0), 1);
z = wb * (wh/wb).^((k + N + 0.5*(1 - gam)) / n);
p = wb * (wh/wb).^((k + N + 0.5*(1 + gam)) / n);
g = wh^gam;
end

function InitConditions(block)
P = block.DialogPrm(1).Data;
Ts = P(10);
[zi, pim, gi] = oustaloup(1 - P(4));   % voie integrale : ordre 1 - lambda
[zd, pdm, gd] = oustaloup(P(5));       % voie derivee   : ordre mu
coef = [zi; pim; exp(-pim*Ts); (1 - exp(-pim*Ts))./pim; ...
        zd; pdm; exp(-pdm*Ts); (1 - exp(-pdm*Ts))./pdm; gi; gd];
block.Dwork(5).Data = coef;
block.Dwork(1).Data = zeros(7, 1);
block.Dwork(2).Data = zeros(7, 1);
block.Dwork(3).Data = [0; 0];          % [xI ; xF]
% Rapport cyclique initial = feedforward a la consigne initiale (comme le noyau).
R = P(15); rL = P(16); rd = P(17); Ron = P(18); Vin = P(19);
v0 = 24;
dff0 = (R + rL + rd)*v0 / (R*Vin - (Ron - rd)*v0);
block.Dwork(4).Data = [dff0; dff0];    % [d applique ; d_pre associe]
end

function Outputs(block)
d = block.Dwork(4).Data;
P = block.DialogPrm(1).Data;
block.OutputPort(1).Data = d(1);
block.OutputPort(2).Data = d(2);
block.OutputPort(3).Data = P(4);
end

function Update(block)
P = block.DialogPrm(1).Data;
kp = P(1); ki = P(2); kd = P(3); wf = P(6); b = P(7); c = P(8); Kb = P(9);
Ts = P(10); Vn = P(11); Dn = P(12); dmin = P(13); dmax = P(14);
R = P(15); rL = P(16); rd = P(17); Ron = P(18); Vin = P(19);

r  = block.InputPort(1).Data;
vo = block.InputPort(2).Data;

C  = block.Dwork(5).Data;
zi = C(1:7);   pim = C(8:14);  phii = C(15:21); betai = C(22:28);
zd = C(29:35); pdm = C(36:42); phid = C(43:49); betad = C(50:56);
gi = C(57); gd = C(58);

wi = block.Dwork(1).Data; wd = block.Dwork(2).Data; xs = block.Dwork(3).Data;
xI = xs(1); xF = xs(2);

en = (r - vo)/Vn; eb = (b*r - vo)/Vn; ec = (c*r - vo)/Vn;

y = en;                                  % voie integrale
for i = 1:7
    yi = y + (zi(i) - pim(i))*wi(i);
    wi(i) = phii(i)*wi(i) + betai(i)*y;
    y = yi;
end
q = gi*y;

phiF = exp(-wf*Ts);                      % voie derivee
xF = phiF*xF + (1 - phiF)*ec;
y = xF;
for i = 1:7
    yi = y + (zd(i) - pdm(i))*wd(i);
    wd(i) = phid(i)*wd(i) + betad(i)*y;
    y = yi;
end
dterm = gd*y;

u = kp*eb + xI + kd*dterm;
dff = (R + rL + rd)*r / (R*Vin - (Ron - rd)*r);
d_tilde = dff + Dn*u;
d_sat = min(max(d_tilde, dmin), dmax);
xI = xI + Ts*(ki*q + Kb*(d_sat - d_tilde)/Dn);

block.Dwork(1).Data = wi; block.Dwork(2).Data = wd;
block.Dwork(3).Data = [xI; xF];
% Retard d'un echantillon : la sortie de l'instant suivant est d_sat calcule ici.
block.Dwork(4).Data = [d_sat; d_tilde];
end
