% First panel of phase-shift section figure :

% For a few exemples with very small frequency difference, show that the trajectory is not
% circular anymore.

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
params.l2 = 2 / 1.5 * flagLength;
assignForceAndTorqueCoeffs;
params.phis = [pi/6;-pi/6];

% For generating various beat patterns.
% coeff = rand(1,4);
% alphaFun = @(t) 2*coeff(1)*sin(2*pi*t*baseFreq) - 2*coeff(2);
% betaFun = @(t) 2*coeff(3)*cos(2*pi*t*baseFreq + 9*pi/10) - 2*coeff(4);

% Chlamy-like beat pattern.
alphaFun = @(t) sin(2*pi*t*baseFreq) - 0.5;
betaFun = @(t) cos(2*pi*t*baseFreq + 9*pi/10) - 1;

params.alphas = @(t) alphaFun(t) * [1;-1];
params.betas = @(t) betaFun(t) * [1;-1];

init = [0;0;0];

figure(4);clf;
tt=0:0.01:1;
plot(alphaFun(tt),betaFun(tt))

% Select a few freq multipliers and sweep through them.

ts = linspace(0,200,1e4);

freqMultipliers = [1.1 1.05 1.02 1.01 1.005 1.002 1.0005];
allPositions = cell(length(freqMultipliers),1);
allPhis = cell(length(freqMultipliers),1);
allForces = cell(length(phaseShifts),1);

for ind = 1 : length(freqMultipliers)
    disp(ind)
    freqMultiplier = freqMultipliers(ind);

    params.alphas = @(t) [alphaFun(freqMultiplier*t); -alphaFun(t)];
    params.betas = @(t) [betaFun(freqMultiplier*t); -betaFun(t)];

    [positions, phis, contact_forces] = solve_problem(params,init,ts);
    % Save the positions.
    allPositions{ind} = positions;
    allPhis{ind} = phis;
    allForces{ind} = contact_forces;

    if videos
        make_video(positions, phis, ts, params, 30, 1, ['PS',num2str(phaseShift)], true);
    end

end

%% Plots.

% A plot with the x and y components of each force.
figure(1);clf;
tiledlayout(min(6,length(freqMultipliers)),2)

Favg = zeros(length(freqMultipliers),2);
for i = 1:length(freqMultipliers)
    F = allForces{i};
    Favg(i,1) = sum(F(1,:)-F(3,:))/length(F(1,:));
    Favg(i,2) = sum(F(2,:)+F(4,:))/length(F(2,:));
end


    for i = 1:min(6,length(freqMultipliers))
        nexttile(2*i-1);
        F = allForces{i};
        plot(ts(1:end-1),F(1,:),'b');
        hold on
        plot(ts(1:end-1),F(3,:),'r');
    
        nexttile(2*i);
        plot(ts(1:end-1),F(2,:),'b');
        hold on
        plot(ts(1:end-1),-F(4,:),'r');
    end


%%

% A figure with the trajectories/
figure(3);clf;
set(gcf,'color','w');
set(gcf, 'Position',  [1,1,600,550])

hold on
box on
colororder(othercolor('Paired6',length(freqMultipliers)));
for ind = 1 : length(freqMultipliers)
    plot(allPositions{ind}(:,1),allPositions{ind}(:,2),'DisplayName',num2str(freqMultipliers(ind)))
end

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
axis equal

leg = legend;
leg.Interpreter = 'latex';
leg.Position = [0.2 0.24 0.1 0.2];
set(gca,'TickLabelInterpreter','latex') 
set(gca,'FontSize',22) 

exportgraphics(gcf,'fig_phase_panel1.pdf','ContentType','vector')

