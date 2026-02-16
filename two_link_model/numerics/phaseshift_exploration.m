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

coeff = rand(1,4);
alphaFun = @(t) 2*coeff(1)*sin(2*pi*t*baseFreq) - 2*coeff(2);
betaFun = @(t) 2*coeff(3)*cos(2*pi*t*baseFreq + 9*pi/10) - 2*coeff(4);

% alphaFun = @(t) sin(2*pi*t*baseFreq) - 0.5;
% betaFun = @(t) cos(2*pi*t*baseFreq + 9*pi/10) - 1;

%alphaFun = @(t) sin(2*pi*t*baseFreq);
%betaFun = @(t) cos(2*pi*t*baseFreq + 9*pi/10);

params.alphas = @(t) alphaFun(t) * [1;-1];
params.betas = @(t) betaFun(t) * [1;-1];

init = [0;0;0];

figure(4);clf;
tt=0:0.01:1;
plot(alphaFun(tt),betaFun(tt))


% Sweep through phase shifts, and compute the average force.

ts = linspace(0,20,1e4);
phaseShifts = 0:0.05:1;
allPositions = cell(length(phaseShifts),1);
allPhis = cell(length(phaseShifts),1);
allForces = cell(length(phaseShifts),1);

for ind = 1 : length(phaseShifts)
    disp(ind)
    phaseShift = phaseShifts(ind);

    params.alphas = @(t) [alphaFun(t); -alphaFun(t+phaseShift)];
    params.betas = @(t) [betaFun(t); -betaFun(t+phaseShift)];

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
tiledlayout(min(6,length(phaseShifts)),2)

Favg = zeros(length(phaseShifts),2);
for i = 1:length(phaseShifts)
    F = allForces{i};
    Favg(i,1) = sum(F(1,:)-F(3,:))/length(F(1,:));
    Favg(i,2) = sum(F(2,:)+F(4,:))/length(F(2,:));
end


    for i = 1:min(6,length(phaseShifts))
        nexttile(2*i-1);
        F = allForces{i};
        % plot(ts(1:end-1),F(1,:),'b');
        hold on
        plot(ts(1:end-1),F(1,:).^2 + F(2,:).^2,'g');
        plot(ts(1:end-1),F(3,:).^2 + F(4,:).^2,'black');
    
        nexttile(2*i);
        % plot(ts(1:end-1),F(2,:),'b');
        hold on
        % plot(ts(1:end-1),-F(4,:),'r');
    end


% Linear regression on the phi dynamics to get the average angular
% velocity.
phidotAvg = zeros(length(phaseShifts),1);
for ind = 1 : length(phaseShifts)
    P = polyfit(ts,allPhis{ind},1);
    phidotAvg(ind) = P(1);
end

% A plot of key values.
figure(2);clf;
tiledlayout(1,3)

% Average difference in x.
nexttile(1);
plot(phaseShifts,Favg(:,1));
grid on

% Average difference in y.
nexttile(2);
plot(phaseShifts,Favg(:,2));
grid on

% Angular velocity. 
nexttile(3);
plot(phaseShifts,phidotAvg,'Color','black','LineWidth',2)
xlabel('Phase Shifts')
ylabel('Average angular speed')
grid on

%% 
figure(4);clf;
% Try a linear combination of averages of Fx and Fy.
Fcomb = -1.8*Favg(:,1)-Favg(:,2);
plot(phaseShifts,Fcomb);
hold on;
PhiM = max(abs(phidotAvg));
FM = max(abs(Fcomb));
plot(phaseShifts,phidotAvg/PhiM*FM)
grid on
axis tight
title(num2str(coeff))

%%

% A figure with the trajectories/
figure(3);clf;
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

