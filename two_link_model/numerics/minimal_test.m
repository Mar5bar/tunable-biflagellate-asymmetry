% theta = @(t) -A*(cos(w*t) - cos(t) + (w-1)*D*t);

[ts,Y] = ode15s(@odefun, [0,100*pi], [0;0;0]);
plot(Y(:,2),Y(:,3))
axis equal

function dY = odefun(t,Y)
    w = 1.5;
    A = 0.1;
    D = 0;
    Dv = 0;
    V = 0.1;
    theta = Y(1);
    dY = zeros(3,1);
    f = @(t) sin(t)+0.1;
    dY(1) = A*(w*f(w*t) - f(t) + (w-1)*D);
    dY(2) = V*(w*f(w*t) + f(t) + (w+1)*Dv)*cos(theta);
    dY(3) = V*(w*f(w*t) + f(t) + (w+1)*Dv)*sin(theta);
end