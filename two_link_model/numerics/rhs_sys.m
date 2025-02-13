function rhs = rhs_sys(alphas, betas, params)
%% Form the linear system that computes the forces and torques due to
%swimmer-frame deformation in the swimmer frame from the params and the
%current % state. Entries are copied directly from symbolic algebra.
    alpha_1 = alphas(1);
    alpha_2 = alphas(2);

    beta_1 = betas(1);
    beta_2 = betas(2);

    phi1 = params.phis(1);
    phi2 = params.phis(2);

    mu = params.mu;
    Cpar = params.Cpar;
    Cperp = params.Cperp;
    r = params.r;
    l1 = params.l1;
    l2 = params.l2;

    rhs = zeros(3,4);

    rhs(1,1) = (Cperp*l1^2*sin(alpha_1 - phi1))/2 + Cpar*l1*l2*(sin(alpha_1 - phi1)/2 - sin(alpha_1 + 2*beta_1 - phi1)/2) + (Cperp*l2*sin(alpha_1 + beta_1 - phi1)*(l2 + 2*l1*cos(beta_1)))/2;
    rhs(1,2) = (Cperp*l1^2*sin(alpha_2 - phi2))/2 + Cpar*l1*l2*(sin(alpha_2 - phi2)/2 - sin(alpha_2 + 2*beta_2 - phi2)/2) + (Cperp*l2*sin(alpha_2 + beta_2 - phi2)*(l2 + 2*l1*cos(beta_2)))/2;
    rhs(1,3) = (Cperp*l2^2*sin(alpha_1 + beta_1 - phi1))/2;
    rhs(1,4) = (Cperp*l2^2*sin(alpha_2 + beta_2 - phi2))/2;
    rhs(2,1) = (Cperp*l1^2*cos(alpha_1 - phi1))/2 - Cpar*l1*l2*(cos(alpha_1 + 2*beta_1 - phi1)/2 - cos(alpha_1 - phi1)/2) + (Cperp*l2*cos(alpha_1 + beta_1 - phi1)*(l2 + 2*l1*cos(beta_1)))/2;
    rhs(2,2) = (Cperp*l1^2*cos(alpha_2 - phi2))/2 - Cpar*l1*l2*(cos(alpha_2 + 2*beta_2 - phi2)/2 - cos(alpha_2 - phi2)/2) + (Cperp*l2*cos(alpha_2 + beta_2 - phi2)*(l2 + 2*l1*cos(beta_2)))/2;
    rhs(2,3) = (Cperp*l2^2*cos(alpha_1 + beta_1 - phi1))/2;
    rhs(2,4) = (Cperp*l2^2*cos(alpha_2 + beta_2 - phi2))/2;
    rhs(3,1) = Cperp*(l2^3/3 + (l1*l2^2*cos(beta_1))/2 + (l2*(r*cos(alpha_1 + beta_1) + l1*cos(beta_1))*(l2 + 2*l1*cos(beta_1)))/2) + (Cperp*l1^2*(2*l1 + 3*r*cos(alpha_1)))/6 + Cpar*l1*l2*sin(beta_1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1));
    rhs(3,2) = Cperp*(l2^3/3 + (l1*l2^2*cos(beta_2))/2 + (l2*(r*cos(alpha_2 + beta_2) + l1*cos(beta_2))*(l2 + 2*l1*cos(beta_2)))/2) + (Cperp*l1^2*(2*l1 + 3*r*cos(alpha_2)))/6 + Cpar*l1*l2*sin(beta_2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2));
    rhs(3,3) = Cperp*((l2^2*(r*cos(alpha_1 + beta_1) + l1*cos(beta_1)))/2 + l2^3/3);
    rhs(3,4) = Cperp*((l2^2*(r*cos(alpha_2 + beta_2) + l1*cos(beta_2)))/2 + l2^3/3);

end