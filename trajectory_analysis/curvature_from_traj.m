function curvature = curvature_from_traj(x,y)

    % Fit a circle.
    par = circleFit(x,y);
    % Curvature is 1/radius.
    curvature = 1/par(3);

end