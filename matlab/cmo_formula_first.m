% CMO_FORMULA_FIRST  Independent MATLAB check of the analytical chain.
%   Eqs. (7)-(12), CARE (20), Bryson weights (21), exact PIDF map (33),
%   coefficient and frequency-domain residuals, and the Proposition 3 identity
%   x_F(t) = c r - v_C(t). Compare with results/cmo/formula_first.json.
Vin = 48; L = 100e-6; C = 100e-6; R = 10; rL = 0.1; rC = 0.05; Ron = 0.05; rd = 0.02; Vf = 0.42;
Ts = 10e-6; dmin = 0.02; dmax = 0.98;

Vo = 24; IL = Vo / R;
d0 = (Vo + (rL + rd) * IL + Vf) / (Vin - (Ron - rd) * IL + Vf);
g = R / (R + rC); rho0 = rL + rd + d0 * (Ron - rd); Ed = Vin + Vf - (Ron - rd) * IL;
A = [-(rho0 + g * rC) / L, -g / L; g / C, -g / (R * C)];
Bd = [Ed / L; 0]; Cy = [g * rC, g];
Gvd = tf(Ed * R * [rC * C 1], [L * C * (R + rC), L + C * ((R + rC) * rho0 + R * rC), R + rho0]);
fprintf('state-space vs closed-form Gvd mismatch: %.3e\n', norm(ss(A, Bd, Cy, 0) - Gvd, inf) / norm(Gvd, inf));

xi0 = log10([15 46 2 * pi * 1500 2 * pi * 1500]);
Di = 10^xi0(1); Dv = 10^xi0(2); wB = 10^xi0(3); Kb = 10^xi0(4); Dd = (dmax - dmin) / 2;
Q = diag([Di^-2, Dv^-2, (Dv / wB)^-2]); Rd = Dd^-2;
Aa = [A zeros(2, 1); -Cy 0]; Ba = [Bd; 0];
[P, K] = icare(Aa, Ba, Q, Rd);
res = Aa' * P + P * Aa - P * Ba / Rd * Ba' * P + Q;
fprintf('K = [%.10g %.10g %.10g], CARE residual %.3e (rel %.3e)\n', K, norm(res), norm(res) / norm(Q));
disp('closed-loop poles:'); disp(eig(Aa - Ba * K));

Kp = K(2) + K(1) / R; Ki = -K(3); Kd = C * (K(1) - rC * K(2)); N = 1 / (rC * C);
theta = [Kp Ki Kd N 1 0 Kb];
fprintf('theta = [%s]\n', num2str(theta, 10));
s = tf('s');
Clqi = minreal(K(1:2) * ((s * eye(2) - A) \ Bd) / (Cy * ((s * eye(2) - A) \ Bd)) - K(3) / s);
Cpidf = Kp + Ki / s + Kd * N * s / (s + N);
w = logspace(-1, 8, 4000);
e = abs(squeeze(freqresp(Clqi - Cpidf, w))) ./ abs(squeeze(freqresp(Clqi, w)));
fprintf('max relative frequency-domain error: %.3e\n', max(e));

% Proposition 3: with xF' = N (c r - vo - xF), eps = xF - (c r - vC) obeys eps' = -N eps
% whatever the inductor branch, duty saturation or conduction mode (ip = 0).
KdN = Kd * N;
fprintf('identity Kix/rC - Kvx = Kd N : %.3e\n', abs((K(1) / rC - K(2)) - KdN));
fprintf('N < pi/Ts : %d  (N = %.4g, pi/Ts = %.4g)\n', N < pi / Ts, N, pi / Ts);
