addpath(genpath('../../common'))

% Make videos?
videos = false;

% We'll run through beats where we alter the frequency of one of the flagella.
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

init = [0;0;0];

% Make a nice video to show what's happening over a few beats.
ts = linspace(0,50,5e4);
[positions, phis] = solve_problem(params,init,ts);
if videos
    make_video(positions, phis, ts, params, 30, 1, 'videos/normal_swimming', true, false);
end

% Simulate for a long time to see differences in messy trajectories.
ts = linspace(0,50,1e4);
freqMultipliers = 1.1:0.1:5;
allPositions = cell(length(freqMultipliers),1);
allPhis = cell(length(freqMultipliers),1);
for ind = 1 : length(freqMultipliers)
    disp(ind)
    freqMultiplier = freqMultipliers(ind);

    params.alphas = @(t) [alphaFun(freqMultiplier*t); -alphaFun(t)];
    params.betas = @(t) [betaFun(freqMultiplier*t); -betaFun(t)];

    [positions, phis] = solve_problem(params,init,ts);
    % Save the positions.
    allPositions{ind} = positions;
    allPhis{ind} = phis;

    if videos
        make_video(positions, phis, ts, params, 30, 1, ['videos/freqMultiplier',num2str(freqMultiplier)], true, false);
    end
end

figure
font_size = 20;
hold on
box on
colororder(viridis(length(freqMultipliers)));
for ind = 1 : length(freqMultipliers)
    plot(allPositions{ind}(:,1),allPositions{ind}(:,2))
end
c = colorbar;
c.Label.String = "Frequency multiplier";
c.Label.Interpreter = "latex";
c.TickLabelInterpreter = "latex";
colormap(viridis)
caxis([min(freqMultipliers), max(freqMultipliers)])
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
axis equal
exportgraphics(gcf,'trajectory_with_frequency_multiplier.pdf')

% Fit circles to get the radius of curvature as a function of frequency.
radii = zeros(length(freqMultipliers),1);
for ind = 1 : length(freqMultipliers)
    par = circleFit(allPositions{ind}(:,1),allPositions{ind}(:,2));
    radii(ind) = par(3);
end
figure
set(gcf, 'Position',  [611,895,500,400])
set(gcf,'color','w');
hold on
box on
plot(freqMultipliers,radii(1)./radii,'+','Color','black','LineWidth',2)
xlabel('$k_R / k_L$','FontSize',font_size,'Interpreter','latex');
ylabel('$\kappa$','FontSize',font_size,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');   
set(gca,'FontSize',font_size);

fitfun = fittype(@(a,x) a * (x-1)./(x+1))
[fitted_curve, gof] = fit(freqMultipliers', radii(1)./radii, fitfun, 'StartPoint', 1);
a = coeffvalues(fitted_curve);

plot(freqMultipliers,a*(freqMultipliers-1)./(freqMultipliers+1),'k:','LineWidth',1)
ylim([0,15])
grid on
legend({'Simulation','Scaling law'},'location','southeast')
exportgraphics(gcf,'curvature_with_frequency_multiplier.pdf')

save("output.mat")
