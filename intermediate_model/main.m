% Solve the intermediate model and plot the trajectory.

asymmetry_type = "amplitude"; alpha = 1.1;
% asymmetry_type = "frequency"; k2 = 1.1;
% asymmetry_type = "phase"; phi = pi;

% Specify the system.
J = @(t) sin(t);
k1 = 1;
A1 = 1;
A2 = 1;
B1 = 1;
B2 = 1;

T = 1e4;
mu = 1;
radius = 1;
psi = pi/6;

% These are applicable for all types of asymmetry.
f11 = @(t) k1*A1 + k1*B1*J(k1*t);
f12 = @(t) k1*A2 + k1*B2*J(k1*t);

% For amplitude asymmetry:
switch asymmetry_type
    case "amplitude"
        f21 = @(t) k1*alpha*A1 + k1*alpha*B1*J(k1*t);
        f22 = @(t) k1*alpha*A2 + k1*alpha*B2*J(k1*t);
    case "frequency"
        f21 = @(t) k2*A1 + k2*B1*J(k2*t);
        f22 = @(t) k2*A2 + k2*B2*J(k2*t);
    case "phase"
        f21 = @(t) k1*A1 + k1*B1*J(k1*t + phi);
        f22 = @(t) k1*A2 + k1*B2*J(k1*t + phi);
    otherwise
        error("Invalid asymmetry type.");
end

params = {};
params.f11 = f11;
params.f12 = f12;
params.f21 = f21;
params.f22 = f22;
params.radius = radius;
params.mu = mu;
params.psi = psi;

% Solve the ODE system.
ts = linspace(0,T,1e3);
z0 = [0;0;0];
[t,z] = ode15s(@(t,z) ode_rhs(t,z,params), ts, z0);
x = z(:,1);
y = z(:,2);
theta = z(:,3);

plot(x,y);
xlabel('$x$');
ylabel('$y$');
axis equal



function dz = ode_rhs(t, z, params)
    f11 = params.f11;
    f12 = params.f12;
    f21 = params.f21;
    f22 = params.f22;

    [e1x, e1y] = e1_from_theta(z(3));
    [e2x, e2y] = e2_from_theta(z(3));
    u = (f11(t) + f21(t))*e1x + (f12(t) - f22(t))*e2x;
    v = (f11(t) + f21(t))*e1y + (f12(t) - f22(t))*e2y;
    omega = (-f11(t) + f21(t))*sin(params.psi) + (f12(t) - f22(t))*cos(params.psi);

    u = u / (6*pi*params.mu*params.radius);
    v = v / (6*pi*params.mu*params.radius);
    omega = omega / (8*pi*params.mu*params.radius^2);

    dz = [u;v;omega];
end

function [x,y] = e1_from_theta(theta)
    x = cos(theta);
    y = sin(theta);
end

function [x,y] = e2_from_theta(theta)
    x = -sin(theta);
    y = cos(theta);
end