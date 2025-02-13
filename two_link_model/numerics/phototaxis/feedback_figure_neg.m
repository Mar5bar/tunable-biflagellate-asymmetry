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

lightPos = [0,0];

params.lightPos = lightPos;

%% Controller.

% Simple principle : if the light faces the sensor, don't activate the
% differential freq; if it does not, do.
% We'll have to solve the ODE on each time period, and concatenate, to
% avoid discontinuities.

% Frequency multiplier
freq_max = 1.5;
params.freq_max = freq_max;

% Critical angle
a_crit = pi - 1.5*pi/10;
params.a_crit = a_crit;

% Gait function.
alphaFun = @(t) baseAmp*(sin(2*pi*t)-0.5);
betaFun = @(t) baseAmp*(cos(2*pi*t + 9*pi/10)-1);

params.alpha = alphaFun;
params.beta = betaFun;

% Final time.
T = 200;


% Loop over positions.
r = 3e-5; 
N_init = 15;
th_init_var = [0:2*pi/N_init:2*pi*(N_init-1)/N_init]+1e-2;

all_tps = cell(length(th_init_var),1);
all_positions = cell(length(th_init_var),1);
all_phis = cell(length(th_init_var),1);

for i = 1:length(th_init_var)

init_current = [r*cos(th_init_var(i));r*sin(th_init_var(i));0];
positions = [init_current(1),init_current(2)];
phis = init_current(3);

% Time counter.
tps = 0;

% Stop criterion on position.
stop_crit = 0;

% Start the solver loop.
tic
while tps(end) < T && stop_crit < 5e-5

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
    
    init_current = [positions(end,:),phis(end)];
    [positionsTemp, phisTemp] = solve_problem(params,init_current,tpsTemp);
    
    % Concatenate.
    tps = [tps,tpsTemp(2:end)];
    positions = [positions; positionsTemp(2:end,:)];
    phis = [phis; phisTemp(2:end)];

    % Compute stopping criterion
    stop_crit = norm(positions(end,:));
end
toc

% Save the data.
all_tps{i} = tps;
all_positions{i} = positions;
all_phis{i} = phis;

end

%%

addpath('../')

col = colormap(othercolor('Paired6',length(th_init_var)));

figure(3);clf;

set(gcf, 'Position',  [450,1,400,400])
set(gcf,'color','w');

hold on
box on

for i = 1:length(th_init_var)
    plot(all_positions{i}(:,1),all_positions{i}(:,2),'Color',col(i,:),'LineWidth',1)
    hold on

    plot(all_positions{i}(1,1),all_positions{i}(1,2),'k.','MarkerSize',20)
end

    plot(lightPos(1),lightPos(2),'ko','MarkerSize',12,'LineWidth',3)

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
axis equal tight

set(gca,'TickLabelInterpreter','latex');   
set(gca,'FontSize',18);

exportgraphics(gcf,'fig_phototaxis_neg.pdf','ContentType','vector')

