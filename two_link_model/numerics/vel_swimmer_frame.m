function [vel, A, rhs] = vel_swimmer_frame(t, params)
%% Compute the derivative of the lab frame position and orientation.

%% Compute the velocities in the swimmer frame from the controls.
% Evaluate the current swimmer frame configuration, which is prescribed.
    alphas = params.alphas(t);
    betas = params.betas(t);

    alphadot = diffByT(params.alphas, t);
    betadot = diffByT(params.betas, t);

    A = lin_sys(alphas, betas, params);
    rhs = rhs_sys(alphas, betas, params) * [alphadot; betadot];

    vel = A \ rhs;

end