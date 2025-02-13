function t = compute_t(i, params)
    t = [sin(params.phis(i));
         -cos(params.phis(i));
         0];
end