% Solve the minimal model and plot the trajectory.

% Define parameters.
params = struct();
A_1 = 1; params.A_1 = A_1; % base amplitude of first flagellum
A_2 = 1.0; params.A_2 = A_2;% differential amplitude of second flagellum
eta_t = 1/A_1; params.eta_t = eta_t;% translation coefficient
eta_r = 1/A_1; params.eta_r = eta_r;% rotation coefficient
k_1 = 1; params.k_1 = k_1; % frequency of first flagellum
k_2 = 1; params.k_2 = k_2; % frequency of second flagellum
f = @(t) 0.9+0.1*cos(1*t); params.f = f;% periodic activation function
g = @(t) f(t)-0.01*sin(t)^2;%7/11+1.4*sin(2.4*t); 
params.g = g;% different activation function for second flagellum
phi = 0; params.phi = phi; % phase shift between both flagella

% ODE conditions setup.
init = [0;0;0];
T = 140; nb_dt = 3e3;
tps = linspace(0,T,nb_dt);

% ODE solver setup.
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
dZ= @(t,z) dynamics(t,z,params);

% Solving the dynamics.
[tps, traj] = ode45(dZ,tps,init,opts);

% Plot
figure(1);clf;
plot(traj(:,1),traj(:,2),'b','LineWidth',3)
axis equal tight
% figure(2);clf;
% tiledlayout(3,1)
% nexttile(1);
% plot(tps,traj(:,1));
% nexttile(2);
% plot(tps,traj(:,2));
% nexttile(3);
% plot(tps,traj(:,3)+(k_2-k_1)*tps);

function zdot = dynamics(t,z,params)

    % Unpack parameters.
    eta_t = params.eta_t;
    eta_r = params.eta_r;
    A_1 = params.A_1;
    A_2 = params.A_2;
    f = params.f;
    g = params.g;
    k_1 = params.k_1;
    k_2 = params.k_2;
    phi = params.phi;

    % Unpack state.
    x = z(1);
    y = z(2);
    theta = z(3);

    % Define velocity.
    V = eta_t * (A_1 * k_1 * f(k_1*t) + A_2 * k_2 * g(k_2*t + phi));

    % Compute the dynamics.
    zdot = zeros(3,1);
    zdot(1) = V*cos(theta);
    zdot(2) = V*sin(theta);
    zdot(3) = eta_r * (A_1 * k_1 * f(k_1*t) - A_2 * k_2 * g(k_2*t + phi));
end
