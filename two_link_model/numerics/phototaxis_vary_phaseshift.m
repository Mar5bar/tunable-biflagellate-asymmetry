% Make videos?
videos = true;

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
params.l2 = 2 / 1.5 * flagLength;
assignForceAndTorqueCoeffs;
params.phis = [pi/6;-pi/6];

alphaFun = @(t) sin(2*pi*t*baseFreq) - 0.5;
betaFun = @(t) cos(2*pi*t*baseFreq + 9*pi/10) - 1;

%alphaFun = @(t) sin(2*pi*t*baseFreq);
%betaFun = @(t) cos(2*pi*t*baseFreq + 9*pi/10);

params.alphas = @(t) alphaFun(t) * [1;-1];
params.betas = @(t) betaFun(t) * [1;-1];

init = [0;0;0];

figure(4);clf;
tt=0:0.01:1;
plot(alphaFun(tt),betaFun(tt))

%%

% Make a nice video to show what's happening over a few beats.
ts = linspace(0,30,1e4);
figure(2);clf;
[positions, phis] = solve_problem(params,init,ts);
if videos
    make_video(positions, phis, ts, params, 30, 1, 'videos/normal_swimming', true);
end

% Simulate for a long time to see differences in messy trajectories.
ts = linspace(0,10,1e4);
phaseShifts = [0.1 0.2 0.5];
allPositions = cell(length(phaseShifts),1);
allPhis = cell(length(phaseShifts),1);
for ind = 1 : length(phaseShifts)
    disp(ind)
    phaseShift = phaseShifts(ind);

    params.alphas = @(t) [alphaFun(t); -alphaFun(t+phaseShift)];
    params.betas = @(t) [betaFun(t); -betaFun(t+phaseShift)];

    [positions, phis] = solve_problem(params,init,ts);
    % Save the positions.
    allPositions{ind} = positions;
    allPhis{ind} = phis;

    if videos
        make_video(positions, phis, ts, params, 30, 1, ['videos/PS',num2str(phaseShift)], true);
    end
end

%%

figure(1);clf;
hold on
box on
colororder(jet(length(phaseShifts)));
for ind = 1 : length(phaseShifts)
    plot(allPositions{ind}(:,1),allPositions{ind}(:,2))
end
c = colorbar;
c.Label.String = "Frequency multiplier";
c.Label.Interpreter = "latex";
c.TickLabelInterpreter = "latex";
colormap("jet")
caxis([min(phaseShifts), max(phaseShifts)])
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
axis equal
%exportgraphics(gcf,'trajectory_with_frequency_multiplier.png','Resolution',400)

% Linear regression on the phi dynamics to get the average angular
% velocity.

% phidotAvg = zeros(length(phaseShifts),1);
% for ind = 1 : length(phaseShifts)
%     P = polyfit(ts,allPhis{ind},1);
%     phidotAvg(ind) = P(1);
% end
% 
% figure(3);clf;
% box on
% %P = polyfit(phaseShifts,phidotAvg,1);
% plot(phaseShifts,phidotAvg,'Color','black','LineWidth',2)
% xlabel('Phase Shifts')
% ylabel('Average angular speed')
% 
% 
% % Fit circles to get the radius of curvature as a function of frequency.
% radii = zeros(length(phaseShifts),1);
% for ind = 1 : length(phaseShifts)
%     par = circleFit(allPositions{ind}(:,1),allPositions{ind}(:,2));
%     radii(ind) = par(3);
% end
% figure(2);clf;
% box on
% plot(phaseShifts,1./radii,'Color','black','LineWidth',2)
% xlabel('Phase Shift')
% ylabel('Curvature (1/m)')
% %exportgraphics(gcf,'curvature_with_frequency_multiplier.png','Resolution',400)
% 
% save("output.mat")
