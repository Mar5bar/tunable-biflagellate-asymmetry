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


% Sweep through phase shifts, and compute the average force.

ts = linspace(0,50,1e4);
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

%%

% A figure with the trajectories/
figure(3);clf;
set(gcf,'color','w');
set(gcf, 'Position',  [1,1,600,550])

hold on
box on
colororder(othercolor('Paired6',length(phaseShifts)));
for ind = 1 : length(phaseShifts)
    plot(allPositions{ind}(:,1),allPositions{ind}(:,2))
end
colormap(othercolor('Paired6',length(phaseShifts)));
c = colorbar;
c.Label.String = "$\psi$";
c.Label.Rotation = 0;
c.Label.Position = [0.8 1.07];
c.Label.Interpreter = "latex";
c.FontSize = 24;
c.Ticks = [0 0.25 0.5 0.75 1];
c.TickLabels = {'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$'};
c.TickLabelInterpreter = "latex";


caxis([min(phaseShifts), max(phaseShifts)])
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex') 
set(gca,'FontSize',22) 

axis equal

exportgraphics(gcf,'fig_phase_panel2.pdf','ContentType','vector')

