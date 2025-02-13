addpath(genpath('../'))

% We'll run through beats where we alter the amplitude of one of the flagella.
flagLength = 1e-5; % Flagellum length of 10 micrometres.
flagRadius = 10e-9; % Flagellum radius of 10 nanometres.
epsilon = flagRadius / flagLength;

baseAmp = 1;

params = struct();
params.mu = 1;
params.Cpar = -2*pi*params.mu/(log(2/epsilon) - 0.5); % Gray and Hancock.
params.Cperp = 2*params.Cpar;
params.r = 5e-6; % Radius in metres.
params.l1 = 0.5/1.5 * flagLength;
params.l2 = 1 / 1.5 * flagLength;
assignForceAndTorqueCoeffs;
params.phis = [pi/6;-pi/6];

% Position of the light source.
theta = 1.1*pi;
lightPos = -2*[cos(theta),sin(theta)]*1e-5;

params.lightPos = lightPos;

%% Controller.

% Simple principle : if the light faces the sensor, don't activate the
% differential freq; if it does not, do.
% We solve the ODE cycle by cycle, and concatenate, to
% avoid discontinuities.

% Frequency multiplier
freq_max = 2.5;
params.freq_max = freq_max;

% Critical angle
a_crit = 9*pi/10;
params.a_crit = a_crit;

% Gait function.

alphaFun = @(t) baseAmp*(sin(2*pi*t)-0.5);
betaFun = @(t) baseAmp*(cos(2*pi*t + 9*pi/10)-1);

params.alpha = alphaFun;
params.beta = betaFun;

init_current = [0;0;0];
positions = [init_current(1),init_current(2)];
phis = init_current(3);

% Final time.
T = 400;

% Time counter.
tps = 0;

% Start the loop.
while tps(end) < T

    % Compute the angle between source and sensor.
    L = lightPos - positions(end,:); % Relative position of the light source.
    L = L/norm(L);
    D = [cos(phis(end)),sin(phis(end))]; % Orientation of the swimmer.
    lambda = acos(D*L')*sign(det([D',L'])); % Angle between both vectors.

    % Trigger the phototaxis.
    if abs(lambda) > a_crit
        freqMultiplier = 1;
    else
        freqMultiplier = freq_max;
    end

    % Controls.
    params.alphas = @(t) [alphaFun(t); -alphaFun(freqMultiplier*t)];
    params.betas = @(t) [betaFun(t); -betaFun(freqMultiplier*t)];

    % One period.
    t_per = 1/freqMultiplier;

    % Time interval.
    tpsTemp = linspace(tps(end),tps(end)+t_per,100);

    % Solve.
    tic
    init_current = [positions(end,:),phis(end)];
    [positionsTemp, phisTemp] = solve_problem(params,init_current,tpsTemp);
    toc

    % Concatenate.
    tps = [tps,tpsTemp(2:end)];
    positions = [positions; positionsTemp(2:end,:)];
    phis = [phis; phisTemp(2:end)];
end

% Compute the angle lambda.
L = lightPos/norm(lightPos);
D = [cos(phis),sin(phis)];

lambda = zeros(1,length(phis));
for i = 1:length(lambda)
    lambda(i) = acos(D(i,:)*L')*sign(det([D(i,:)',L']));
end

freqs = ones(1,length(lambda)) + (freq_max-1)*(abs(lambda)>a_crit);

%%

figure(3);clf;
plot(tps,lambda)

alphas = zeros(length(phis),2);
betas = zeros(length(phis),2);

for i = 1:length(phis)
    alphas(i,1) = alphaFun(tps(i));
    alphas(i,2) = - alphaFun((1+freqs(i))*tps(i));
    betas(i,1) = betaFun(tps(i));
    betas(i,2) = - betaFun((1+freqs(i))*tps(i));
end

params.alphas = alphas;
params.betas = betas;

figure(2);clf;
hold on
box on

    plot(positions(:,1),positions(:,2))
    hold on
    plot(lightPos(1),lightPos(2),'o')

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
axis equal

% Video

figure(5);clf;

make_video(positions, phis, tps, params, 80, 1, 'videos/test', true);
