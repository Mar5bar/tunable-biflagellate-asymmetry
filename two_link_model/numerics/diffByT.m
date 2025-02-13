function d = diffByT(fun, t)
    dt = 1e-6;
    d = (fun(t + dt) - fun(t-dt)) / (2*dt);
end