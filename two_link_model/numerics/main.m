flagLength = 1e-5; % Flagellum length of 10 micrometres.
flagRadius = 10e-9; % Flagellum radius of 10 nanometres.
epsilon = flagRadius / flagLength;

baseFreq = 1;

params = struct();
params.mu = 1;
params.Cpar = -2*pi*params.mu/(log(2/epsilon) - 0.5); % Gray and Hancock.
params.Cperp = 2*params.Cpar;
params.r = 5e-6; % Radius in metres.
params.l1 = 0.5/1.5 * flagLength;
params.l2 = 1 / 1.5 * flagLength;
assignForceAndTorqueCoeffs;
params.phis = [pi/6;-pi/6];

alphaFun = @(t) sin(2*pi*t*baseFreq) - 0.5;
betaFun = @(t) cos(2*pi*t*baseFreq + 9*pi/10) - 1;

params.alphas = @(t) alphaFun(t) * [1;-1];
params.betas = @(t) betaFun(t) * [1;-1];

initPos = [0;0];
initPhi = 0;
init = [initPos; initPhi];

tEnd = 5;
ts = linspace(0,tEnd,5001);

[positions, phis] = solve_problem(params,init,ts);

make_video(positions, phis, ts, params, 10, 1, './ani', false, false)