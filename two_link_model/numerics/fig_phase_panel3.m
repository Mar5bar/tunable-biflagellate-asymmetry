% Make videos?
videos = false;

% We'll display the phase shift effect for different beat patterns

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

% Different gait designs.

A1s = [1 2   0.5 2   1];
B1s = [1 0.5 2   2   1];
A2s = [1 0.5 2   0.5 2  ];
B2s = [1 2   0.5 0.5 2  ];
P2s = [9*pi/10 7*pi/10 9*pi/10 13*pi/10 13*pi/10];

phidotAvg = zeros(length(phaseShifts),length(A1s));

% Loop on gait designs.

for j = 1:length(A1s) 

    % Define the gait.
    alphaFun = @(t) A1s(j) * sin(2*pi*t*baseFreq) - A2s(j) * 1;
    betaFun = @(t) B1s(j) * cos(2*pi*t*baseFreq + P2s(j)) - B2s(j) * 1;

    params.alphas = @(t) alphaFun(t) * [1;-1];
    params.betas = @(t) betaFun(t) * [1;-1];

    init = [0;0;0];

    figure(4);clf;
    tt=0:0.01:1;
    plot(alphaFun(tt),betaFun(tt))

% Sweep through phase shifts, and compute the average force.

    ts = linspace(0,20,1e4);
    phaseShifts = 0:0.01:1;
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

    % Linear regression on the phi dynamics to get the average angular
    % velocity.
    
    for ind = 1 : length(phaseShifts)
        P = polyfit(ts,allPhis{ind},1);
        phidotAvg(ind,j) = P(1);
    end

end

%% Plots.

% A plot of key values.
figure(2);clf;
set(gcf,'color','w');
set(gcf, 'Position',  [1,1,600,550])
hold on
box on

% Angular velocity. 

colororder(othercolor('Paired6',length(A1s)));

str_p = {'9\pi/10','7\pi/10','9\pi/10','-7\pi/10','-7\pi/10'};

for j = 1:length(A1s)
    plot(phaseShifts,phidotAvg(:,j),'LineWidth',3,...
        'DisplayName',['$(',num2str(A1s(j)),',',num2str(A2s(j)),',',num2str(B1s(j)),',',num2str(B2s(j)),',',str_p{j},')$'])
    xlabel('Phase Shifts')
    ylabel('Average angular speed')
    grid on
end

xlabel('$\psi$','Interpreter','latex')
ylabel('$\Omega$','Interpreter','latex')

leg = legend;
leg.Interpreter = 'latex';
leg.Location = "northwest";

set(gca,'TickLabelInterpreter','latex') 
set(gca,'FontSize',22) 

legend

exportgraphics(gcf,'fig_phase_panel3.pdf','ContentType','vector')
