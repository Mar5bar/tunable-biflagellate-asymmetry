function [positions, phis] = solve_problem(params,init,ts)
    if nargin < 1
        params = struct();
        params.mu = 1;
        params.Cpar = -1;
        params.Cperp = -2;
        params.eta = 10;
        params.r = 1;
        params.l1 = 1;
        params.l2 = 1;
        assignForceAndTorqueCoeffs;
        params.alphas = @(t) [1*sin(t);
                              -1*sin(t)];
        params.betas = @(t)  [1*cos(t);
                              -1*cos(t)];
        params.phis = [pi/4;-pi/4];
    end
    if nargin < 2
        init = [0;0;0]
    end
    if nargin < 3
        ts = linspace(0,2*pi,1e2);
    end

    opts = odeset('AbsTol',1e-12,'RelTol',1e-9);
    f = @(t,y) ode_fun(t,y,params);
    [ts, sols] = ode15s(f,ts,init,opts);
    positions = sols(:,1:2);
    phis = sols(:,3);
end

function dy = ode_fun(t,y,params)
    phi = y(3);
    velSwimmerFrame = vel_swimmer_frame(t, params);
    velLab = rot_mat(phi) * [velSwimmerFrame(1:2);0];
    dy = [velLab(1:2); velSwimmerFrame(3)];
end 