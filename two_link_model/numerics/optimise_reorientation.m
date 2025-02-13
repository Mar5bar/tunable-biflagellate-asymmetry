% Given a gait, numerically compute the optimal frequency difference that results in the largest per-cycle reorientation of the swimmer.
% For each frequency, we attempt to simulate an integer number of
% beats of both flagella, and the change in angle by the number of beats of the base flagellum.

% We'll run through beats where we alter the frequency of one of the flagella.
flagLength = 1e-5; % Flagellum length of 10 micrometres.
flagRadius = 10e-9; % Flagellum radius of 10 nanometres.
epsilon = flagRadius / flagLength;


params = struct();
params.baseFreq = 1;
params.mu = 1;
params.Cpar = -2*pi*params.mu/(log(2/epsilon) - 0.5); % Gray and Hancock.
params.Cperp = 2*params.Cpar;
params.r = 5e-6; % Radius in metres.
params.l1 = 0.5/1.5 * flagLength;
params.l2 = 1 / 1.5 * flagLength;
assignForceAndTorqueCoeffs;
params.phis = [pi/6;-pi/6];

options = optimset('Display','iter');
[fOpt, reorientationOpt] = fminbnd(@(f) obj_fun(f, params), 1, 10, options);

function reorientation = obj_fun(freqMultiplier, params)
    % Specify the beat here. This can easily be converted to a function
    % that accepts alphaFun and betaFun as arguments.
    alphaFun = @(t) sin(2*pi*t*params.baseFreq) - 0.5;
    betaFun = @(t) cos(2*pi*t*params.baseFreq + 9*pi/10) - 1;

    % Given a frequency multiplier, simulate the swimming and output the reorientation rate.
    params.alphas = @(t) [alphaFun(t); -alphaFun(freqMultiplier*t)];
    params.betas = @(t) [betaFun(t); -betaFun(freqMultiplier*t)];

    % Find when both beats are aligned most closely in phase.
    % Find a rational approximation to freqMultiplier with tolerance 1e-6.
    % The denominator will give us the number of periods of the base frequency
    % to simulate.
    [~,tMax] = rat(freqMultiplier, 1e-2);
    ts = [0, tMax];
    [positions, phis] = solve_problem(params,[0;0;0],ts);
    reorientation = -abs(phis(end) - phis(1)) / tMax;
end