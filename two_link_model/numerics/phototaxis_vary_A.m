addpath(genpath('../../common'))

% Make videos?
videos = false;

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

alphaFun = @(t) baseAmp*(sin(2*pi*t) - 0.5);
betaFun = @(t) baseAmp*(cos(2*pi*t + 9*pi/10) - 1);

params.alphas = @(t) alphaFun(t) * [1;-1];
params.betas = @(t) betaFun(t) * [1;-1];

init = [0;0;0];

% Make a nice video to show what's happening over a few beats.
ts = linspace(0,50,1e4);
[positions, phis] = solve_problem(params,init,ts);
if videos
    make_video(positions, phis, ts, params, 30, 1, 'videos/normal_swimming', true);
end

% Simulate for a long time to see differences in messy trajectories.
ts = linspace(0,100,1e4);
ampMultipliers = linspace(0.5,1.2,21);
allPositions = cell(length(ampMultipliers),1);
allPhis = cell(length(ampMultipliers),1);
for ind = 1 : length(ampMultipliers)
    disp(ind)
    ampMultiplier = ampMultipliers(ind);

    params.alphas = @(t) [alphaFun(t); -ampMultiplier*alphaFun(t)];
    params.betas = @(t) [betaFun(t); -ampMultiplier*betaFun(t)];

    [positions, phis] = solve_problem(params,init,ts);
    % Save the positions.
    allPositions{ind} = positions;
    allPhis{ind} = phis;

    if videos
        make_video(positions, phis, ts, params, 30, 1, ['videos/ampMultiplier',num2str(ampMultiplier)], true);
    end
end

figure
hold on
box on
colororder(viridis(length(ampMultipliers)));
for ind = 1 : length(ampMultipliers)
    plot(allPositions{ind}(:,1),allPositions{ind}(:,2))
end
c = colorbar;
c.Label.String = "amplitude multiplier";
c.Label.Interpreter = "latex";
c.TickLabelInterpreter = "latex";
colormap(viridis)
caxis([min(ampMultipliers), max(ampMultipliers)])
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
axis equal
exportgraphics(gcf,'trajectory_with_amplitude_multiplier.png','Resolution',400)

% Fit circles to get the radius of curvature as a function of amplitude.
radii = zeros(length(ampMultipliers),1);
for ind = 1 : length(ampMultipliers)
    par = circleFit(allPositions{ind}(:,1),allPositions{ind}(:,2));
    radii(ind) = par(3);
end
figure
box on
plot(ampMultipliers,1./radii,'Color','black','LineWidth',2)
xlabel('Amplitude (angular) multiplier')
ylabel('Curvature (1/m)')
exportgraphics(gcf,'curvature_with_amplitude_multiplier.png','Resolution',400)

save("output.mat")
