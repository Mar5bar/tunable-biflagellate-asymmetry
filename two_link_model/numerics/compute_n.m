function n = compute_n(i, params)
    n = [cos(params.phis(i));
         sin(params.phis(i));
         0];
end