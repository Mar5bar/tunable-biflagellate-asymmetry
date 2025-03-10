addpath(genpath('../common'))

% Figure for curvature model. 

% Define parameters.
params = struct();
A_1 = 1; params.A_1 = A_1; % base amplitude of first flagellum
A_2 = 1; params.A_2 = A_2;% differential amplitude of second flagellum
eta_t = 1; params.eta_t = eta_t;% translation coefficient
eta_r = 1; params.eta_r = eta_r;% rotation coefficient
k_1 = 1; params.k_1 = k_1; % frequency of first flagellum
k_2 = 1.3; params.k_2 = k_2; % frequency of second flagellum

% ODE conditions setup.
init = [0;0;0];
T = 25; nb_dt = 1e4;
tps = linspace(0,T,nb_dt);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Different coupling force functions. 
f_var = cell([]);
f_0 = 1;
f_var{1} = @(t) f_0+0.02*sin(t);
f_var{2} = @(t) f_0+0.05*sin(t);
f_var{3} = @(t) f_0+0.1*sin(t);
f_var{4} = @(t) f_0+0.1*sin(t)^2;
f_var{5} = @(t) f_0+0.05*sin(2*t)+0.05*cos(3*t)+0.05*sin(4*t);

% Variable frequency.
k_var = 1.1:0.15:5;

all_tps = cell(length(f_var),length(k_var));
all_x = cell(length(f_var),length(k_var));
all_y = cell(length(f_var),length(k_var));
all_th = cell(length(f_var),length(k_var));

radii = zeros(length(f_var),length(k_var));

% Loop over force functions.

for i_f = 1:length(f_var)

    f = f_var{i_f};params.f = f;

for i = 1:length(k_var)

    k_2 = k_var(i);params.k_2 = k_2;

    % ODE solver setup.
    dZ= @(t,z) dynamics(t,z,params);

    % Solving the dynamics.
    [tps, traj] = ode45(dZ,tps,init,opts);

    % Get average curvature radius.
        par = circleFit(traj(:,1),traj(:,2));
        radii(i_f,i) = par(3);

end

end

%% Plot.

font_size = 20;
figure(2);clf;
set(gcf, 'Position',  [611,895,500,400])
set(gcf,'color','w');
hold on

col = colormap(othercolor('Paired6',length(f_var)));

plot(k_var,((k_var(1)+1)/(k_var(1)-1))*(k_var-1)./(k_var+1),'k:','LineWidth',1)

for i = 1:length(f_var)
    plot(k_var,radii(i,1)./radii(i,:),'+','Color',col(i,:),'LineWidth',2)
end

% axis equal tight
grid on
box on
    xlabel('$k_1 / k_2$','FontSize',font_size,'Interpreter','latex');
    ylabel('$\kappa$','FontSize',font_size,'Interpreter','latex');
    % title('$f(t) = 1 + 0.01 \sin(t)$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');   
    set(gca,'FontSize',font_size);

     l = legend('$(k_1-k_2)/(k_1+k_2)$',...
                '$f(t) = 1+0.02 \sin t$',...
                '$f(t) = 1+0.05 \sin t$',...
                '$f(t) = 1+0.1 \sin t$',...
                '$f(t) = 1+0.1 \sin^2 t$',...
                '$f(t) = g(t) $');
     l.Location = "southeast";
     l.Interpreter = "latex";
     l.FontSize = 18;

exportgraphics(gcf,'fig_phototaxis_fit.pdf','ContentType','vector')

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
