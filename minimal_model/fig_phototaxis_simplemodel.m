% Define parameters.
params = struct();
A_1 = 1; params.A_1 = A_1; % base amplitude of first flagellum
A_2 = 1; params.A_2 = A_2;% differential amplitude of second flagellum
eta_t = 1; params.eta_t = eta_t;% translation coefficient
eta_r = 1; params.eta_r = eta_r;% rotatoon coefficient
k_1 = 1; params.k_1 = k_1; % frequency of first flagellum
k_2 = 1.3; params.k_2 = k_2; % frequency of second flagellum

% ODE conditions setup.
init = [0;0;0];
T = 50; nb_dt = 1e4;
tps = linspace(0,T,nb_dt);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% First activation, simple sine.
f = @(t) 1+0.02*sin(t); params.f = f;% periodic activation function
% Variable frequency.
k_var = [1.5 1.75 1.9 2 2.1 2.5 3 4];

all_tps = cell(2,length(k_var));
all_x = cell(2,length(k_var));
all_y = cell(2,length(k_var));
all_th = cell(2,length(k_var));

for i = 1:length(k_var)
    k_2 = k_var(i);params.k_2 = k_2;
    % ODE solver setup.
    dZ= @(t,z) dynamics(t,z,params);

    % Solving the dynamics.
    [tps, traj] = ode45(dZ,tps,init,opts);

    all_x{1,i} = traj(:,1);
    all_y{1,i} = traj(:,2);
    all_th{1,i} = traj(:,3);
end

% Second activation, more complex function.
f = @(t) 1+0.05*sin(t)^2+0.1*sin(2.5*t); params.f = f;% periodic activation function

for i = 1:length(k_var)
    k_2 = k_var(i);params.k_2 = k_2;
    % ODE solver setup.
    dZ= @(t,z) dynamics(t,z,params);

    % Solving the dynamics.
    [tps, traj] = ode45(dZ,tps,init,opts);

    all_x{2,i} = traj(:,1);
    all_y{2,i} = traj(:,2);
    all_th{2,i} = traj(:,3);
end

% Plot
col = colormap(othercolor('Paired6',length(k_var)));

font_size = 20;
figure(1);clf;
set(gcf, 'Position',  [611,895,500,400])
set(gcf,'color','w');
tl = tiledlayout(2,1);
nexttile(1);
hold on
dx = 5;
for i = 1:length(k_var)
    plot(all_x{1,i}+(i-1)*dx,all_y{1,i},'Color',col(i,:),'LineWidth',1)
end
axis equal tight
grid on
box on
%xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex');
    ylabel('$y$','FontSize',font_size,'Interpreter','latex');
    title('$f(t) = 1 + 0.01 \sin(t)$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');   
    set(gca,'FontSize',font_size);
    c = colorbar;
    c.Position = [0.92 0.185 0.02 0.69];
    c.Ticks = [0:7]/8+1/16;
    c.TickLabels = k_var; 
    c.TickLength = 0;
    c.TickLabelInterpreter='latex';

nexttile(2);
hold on
dx = 5;
for i = 1:length(k_var)
    plot(all_x{2,i}+(i-1)*dx,all_y{2,i},'Color',col(i,:),'LineWidth',1)
end
axis equal tight
ylim([-10 1.5])
grid on
box on
    xlabel('$x$','FontSize',font_size,'Interpreter','latex');
    ylabel('$y$','FontSize',font_size,'Interpreter','latex');
    title('$f(t) = 1 + 0.05 \sin^2(t) + 0.1 \sin (2t)$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');   
    set(gca,'FontSize',font_size);
    yticks([-10 -5 0])

%    tl.Padding="tight";
tl.TileSpacing="tight";

exportgraphics(gcf,'fig_phototaxis_example.pdf','ContentType','vector')

function zdot = dynamics(t,z,params)

    % Unpack parameters.
    eta_t = params.eta_t;
    eta_r = params.eta_r;
    A_1 = params.A_1;
    A_2 = params.A_2;
    f = params.f;
    k_1 = params.k_1;
    k_2 = params.k_2;

    % Unpack state.
    x = z(1);
    y = z(2);
    theta = z(3);

    % Define velocity.
    V = eta_t * (A_1 * k_1 * f(k_1*t) + A_2 * k_2 * f(k_2*t));

    % Compute the dynamics.
    zdot = zeros(3,1);
    zdot(1) = V*cos(theta);
    zdot(2) = V*sin(theta);
    zdot(3) = eta_r * (A_1 * k_1 * f(k_1*t) - A_2 * k_2 * f(k_2*t));
end
