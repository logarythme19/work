function cmo_sfun(block)
% CMO_SFUN  Common sampled controller layer, line-by-line port of cmo/ctrl.py.
%
% Drop-in replacement for buck_comparator_sfun in the block
%   <model>/FO_LQI_fixed_reference Command/Common sampled controller
% with the SAME interface:
%   inputs : 1 Vref, 2 V_o, 3 I_L, 4 I_c, 5 V_in
%   outputs: 1 d_sat (to the PWM, ONE-SAMPLE delay), 2 d_pre, 3 alpha (= 0)
% Single dialog parameter: the 36-element vector p of cmo/ctrl.py
%   p(1)      code: 1 PI, 2 PID, 3 PIDF 2-DOF (CMO-LQI-PIDF), 4 LQR, 5 LQG, 6 LQI
%   p(2:8)    law coefficients
%   p(9:21)   LQR/LQG/LQI matrices
%   p(25:33)  nominal R rL rd Ron Vin Vf rC C Ts  (feedforward, never the plant)
%   p(34:35)  dmin dmax
%   p(36)     initial reference level (0 -> 24 V)
% (MATLAB index = Python index + 1.)

setup(block);
end

function setup(block)
block.NumInputPorts  = 5;
block.NumOutputPorts = 3;
for k = 1:5
    block.InputPort(k).Dimensions = 1;
    block.InputPort(k).DatatypeID = 0;
    block.InputPort(k).Complexity = 'Real';
    block.InputPort(k).SamplingMode = 'Sample';
    block.InputPort(k).DirectFeedthrough = false;
end
for k = 1:3
    block.OutputPort(k).Dimensions = 1;
    block.OutputPort(k).DatatypeID = 0;
    block.OutputPort(k).Complexity = 'Real';
    block.OutputPort(k).SamplingMode = 'Sample';
end
block.NumDialogPrms = 1;
P = block.DialogPrm(1).Data;
block.SampleTimes = [P(33) 0];
block.SimStateCompliance = 'DefaultSimState';
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortSamplingMode);
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitConditions);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Update', @Update);
end

function SetInpPortSamplingMode(block, idx, fd)
block.InputPort(idx).SamplingMode = fd;
for k = 1:block.NumOutputPorts
    block.OutputPort(k).SamplingMode = fd;
end
end

function DoPostPropSetup(block)
block.NumDworks = 1;
block.Dwork(1).Name = 's';
block.Dwork(1).Dimensions = 8;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = 'Real';
block.Dwork(1).UsedAsDiscState = true;
end

function d = ff(r, p)
R = p(25); rL = p(26); rd = p(27); Ron = p(28); Vin = p(29); Vf = p(30);
IL = r / R;
d = (r + (rL + rd) * IL + Vf) / (Vin - (Ron - rd) * IL + Vf);
end

function InitConditions(block)
% Same initialization as ctrl_init: equilibrium at the initial level r0,
% iL(0) = r0/R, vC(0) = r0 (set on the L and C blocks by the runner).
% r0 = p(36) when it is set (cmo_lib base_vector), otherwise 24 V.
p = block.DialogPrm(1).Data;
s = zeros(8, 1);
r0 = 24;
if numel(p) >= 36 && p(36) > 0, r0 = p(36); end
vC0 = r0; iL0 = r0 / p(25); vo0 = r0;
d0 = ff(r0, p);
s(5) = d0; s(6) = d0;
switch p(1)
    case 3, s(2) = p(7) * r0 - vC0;
    case 2, s(2) = -vo0;
    case 5, s(3) = iL0; s(4) = vC0;
end
block.Dwork(1).Data = s;
end

function Outputs(block)
s = block.Dwork(1).Data;
block.OutputPort(1).Data = s(5);
block.OutputPort(2).Data = s(6);
block.OutputPort(3).Data = 0;
end

function Update(block)
p = block.DialogPrm(1).Data;
s = block.Dwork(1).Data;
r = block.InputPort(1).Data; vo = block.InputPort(2).Data;
iL = block.InputPort(3).Data; iC = block.InputPort(4).Data;
code = p(1); Ts = p(33); rC = p(31); R = p(25);
dmin = p(34); dmax = p(35);
dff = ff(r, p);
e = r - vo;
dnow = s(5);
u = 0;
switch code
    case 1
        u = p(2) * e + s(1);
    case 2
        yd = -vo;
        u = p(2) * (p(5) * r - vo) + s(1) + p(4) * (yd - s(2)) / Ts;
        s(2) = yd;
    case 3
        % triangle-hold (FOH) update of the v_C observer x_F with the current
        % sample, then use (Proposition 2: x_F = c r - v_C)
        Kd = p(4); N = p(5); b = p(6); c = p(7);
        w = c * r - vo;
        if s(8) == 0
            s(8) = 1;
        else
            phi = exp(-N * Ts); gam = 1 - (1 - phi) / (N * Ts);
            s(2) = phi * s(2) + (1 - phi) * s(7) + gam * (w - s(7));
        end
        s(7) = w;
        u = p(2) * (b * r - vo) + s(1) + Kd * N * (w - s(2));
    case {4, 6}
        vC = vo - rC * iC;
        u = -p(9) * (iL - r / R) - p(10) * (vC - r);
        if code == 6, u = u - p(11) * s(1); end
    case 5
        xh0 = s(3); xh1 = s(4); xe0 = r / R; xe1 = r;
        yh = r + p(19) * (xh0 - xe0) + p(20) * (xh1 - xe1);
        nu = vo - yh;
        u = -p(9) * (xh0 - xe0) - p(10) * (xh1 - xe1);
        a0 = xh0 + p(11) * nu - xe0; a1 = xh1 + p(12) * nu - xe1;
        s(3) = xe0 + p(13) * a0 + p(14) * a1 + p(17) * (dnow - dff);
        s(4) = xe1 + p(15) * a0 + p(16) * a1 + p(18) * (dnow - dff);
end
dt = dff + u;
ds = min(max(dt, dmin), dmax);
switch code
    case 1, s(1) = s(1) + Ts * (p(3) * e + p(4) * (ds - dt));
    case 2, s(1) = s(1) + Ts * (p(3) * e + p(6) * (ds - dt));
    case 3, s(1) = s(1) + Ts * (p(3) * e + p(8) * (ds - dt));
    case 6, s(1) = s(1) + Ts * (e + p(2) / abs(p(11)) * (ds - dt));
end
s(5) = ds; s(6) = dt;
block.Dwork(1).Data = s;
end
