for i = 1 : 2
    expr = subs(expr, nx(i), cos(phi(i)));
    expr = subs(expr, ny(i), sin(phi(i)));
    expr = subs(expr, nz(i), 0);
    expr = subs(expr, tx(i), sin(phi(i)));
    expr = subs(expr, ty(i), -cos(phi(i)));
    expr = subs(expr, tz(i), 0);
    expr = subs(expr, bx(i), 0);
    expr = subs(expr, by(i), 0);
    expr = subs(expr, bz(i), 1);
end