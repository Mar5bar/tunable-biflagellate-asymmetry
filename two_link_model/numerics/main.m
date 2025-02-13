params = struct();
params.mu = 1;
params.Cpar = -1;
params.Cperp = -2;
params.r = 1;
params.l1 = 1;
params.l2 = 1;

assignForceAndTorqueCoeffs;

alphaFun = @(t) sin(2*pi*t*1) - 0.5;
betaFun = @(t) cos(2*pi*t*1 + 9*pi/10) - 1;

params.alphas = @(t) [alphaFun(2*t);-alphaFun(t)]; % First entry is the "left" flagellum.
params.betas = @(t) [betaFun(2*t);-betaFun(t)];

params.phis = [pi/2;-pi/2];

initPos = [0;0];
initPhi = 0;
init = [initPos; initPhi];

tEnd = 1;
ts = linspace(0,tEnd,1001);

[positions, phis] = solve_problem(params,init,ts);

% make_video(positions, phis, ts, params, 10, 1)