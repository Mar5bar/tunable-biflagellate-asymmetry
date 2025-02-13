function A = lin_sys(alphas, betas, params)
%% Form the linear system that computes the drag on the swimmer in the swimmer
%frame from the params and the current state. Entries are copied directly from
%symbolic algebra.
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

    forceCoeff = params.forceCoeff;
    torqueCoeff = params.torqueCoeff;

    A = zeros(3,3);

    A(1,1) = Cpar*l1 + Cpar*l2 + Cperp*l1 + Cperp*l2 + (Cpar*l2*cos(2*alpha_1 + 2*beta_1 - 2*phi1))/2 + (Cpar*l2*cos(2*alpha_2 + 2*beta_2 - 2*phi2))/2 - (Cperp*l2*cos(2*alpha_1 + 2*beta_1 - 2*phi1))/2 - (Cperp*l2*cos(2*alpha_2 + 2*beta_2 - 2*phi2))/2 + (Cpar*l1*cos(2*alpha_1 - 2*phi1))/2 + (Cpar*l1*cos(2*alpha_2 - 2*phi2))/2 - (Cperp*l1*cos(2*alpha_1 - 2*phi1))/2 - (Cperp*l1*cos(2*alpha_2 - 2*phi2))/2 - forceCoeff*mu*r;
    A(1,2) = -((Cpar - Cperp)*(l1*sin(2*alpha_1 - 2*phi1) + l1*sin(2*alpha_2 - 2*phi2) + l2*sin(2*alpha_1 + 2*beta_1 - 2*phi1) + l2*sin(2*alpha_2 + 2*beta_2 - 2*phi2)))/2;
    A(1,3) = l2*(Cperp*sin(alpha_1 + beta_1 - phi1)*(l2/2 + r*cos(alpha_1 + beta_1) + l1*cos(beta_1)) - Cpar*cos(alpha_1 + beta_1 - phi1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1))) + l2*(Cperp*sin(alpha_2 + beta_2 - phi2)*(l2/2 + r*cos(alpha_2 + beta_2) + l1*cos(beta_2)) - Cpar*cos(alpha_2 + beta_2 - phi2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2))) + l1*((Cperp*l1*sin(alpha_1 - phi1))/2 - Cpar*r*sin(alpha_1)*cos(alpha_1 - phi1) + Cperp*r*sin(alpha_1 - phi1)*cos(alpha_1)) + l1*((Cperp*l1*sin(alpha_2 - phi2))/2 - Cpar*r*sin(alpha_2)*cos(alpha_2 - phi2) + Cperp*r*sin(alpha_2 - phi2)*cos(alpha_2));
    A(2,1) = -((Cpar - Cperp)*(l1*sin(2*alpha_1 - 2*phi1) + l1*sin(2*alpha_2 - 2*phi2) + l2*sin(2*alpha_1 + 2*beta_1 - 2*phi1) + l2*sin(2*alpha_2 + 2*beta_2 - 2*phi2)))/2;
    A(2,2) = Cpar*l1 + Cpar*l2 + Cperp*l1 + Cperp*l2 - (Cpar*l2*cos(2*alpha_1 + 2*beta_1 - 2*phi1))/2 - (Cpar*l2*cos(2*alpha_2 + 2*beta_2 - 2*phi2))/2 + (Cperp*l2*cos(2*alpha_1 + 2*beta_1 - 2*phi1))/2 + (Cperp*l2*cos(2*alpha_2 + 2*beta_2 - 2*phi2))/2 - (Cpar*l1*cos(2*alpha_1 - 2*phi1))/2 - (Cpar*l1*cos(2*alpha_2 - 2*phi2))/2 + (Cperp*l1*cos(2*alpha_1 - 2*phi1))/2 + (Cperp*l1*cos(2*alpha_2 - 2*phi2))/2 - forceCoeff*mu*r;
    A(2,3) = l2*(Cperp*cos(alpha_1 + beta_1 - phi1)*(l2/2 + r*cos(alpha_1 + beta_1) + l1*cos(beta_1)) + Cpar*sin(alpha_1 + beta_1 - phi1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1))) + l2*(Cperp*cos(alpha_2 + beta_2 - phi2)*(l2/2 + r*cos(alpha_2 + beta_2) + l1*cos(beta_2)) + Cpar*sin(alpha_2 + beta_2 - phi2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2))) + l1*((Cperp*l1*cos(alpha_1 - phi1))/2 + Cperp*r*cos(alpha_1)*cos(alpha_1 - phi1) + Cpar*r*sin(alpha_1 - phi1)*sin(alpha_1)) + l1*((Cperp*l1*cos(alpha_2 - phi2))/2 + Cperp*r*cos(alpha_2)*cos(alpha_2 - phi2) + Cpar*r*sin(alpha_2 - phi2)*sin(alpha_2));
    A(3,1) = Cperp*l2*sin(alpha_1 + beta_1 - phi1)*(l2/2 + r*cos(alpha_1 + beta_1) + l1*cos(beta_1)) + Cperp*l2*sin(alpha_2 + beta_2 - phi2)*(l2/2 + r*cos(alpha_2 + beta_2) + l1*cos(beta_2)) + Cpar*l2*cos(alpha_1 + beta_1 - phi1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1)) + Cpar*l2*cos(alpha_2 + beta_2 - phi2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2)) + Cperp*l1*sin(alpha_1 - phi1)*(l1/2 + r*cos(alpha_1)) + Cperp*l1*sin(alpha_2 - phi2)*(l1/2 + r*cos(alpha_2)) - Cpar*l1*r*sin(alpha_1)*cos(alpha_1 - phi1) - Cpar*l1*r*sin(alpha_2)*cos(alpha_2 - phi2);
    A(3,2) = Cperp*l2*cos(alpha_1 + beta_1 - phi1)*(l2/2 + r*cos(alpha_1 + beta_1) + l1*cos(beta_1)) + Cperp*l2*cos(alpha_2 + beta_2 - phi2)*(l2/2 + r*cos(alpha_2 + beta_2) + l1*cos(beta_2)) + Cperp*l1*cos(alpha_1 - phi1)*(l1/2 + r*cos(alpha_1)) + Cperp*l1*cos(alpha_2 - phi2)*(l1/2 + r*cos(alpha_2)) - Cpar*l2*sin(alpha_1 + beta_1 - phi1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1)) - Cpar*l2*sin(alpha_2 + beta_2 - phi2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2)) + Cpar*l1*r*sin(alpha_1 - phi1)*sin(alpha_1) + Cpar*l1*r*sin(alpha_2 - phi2)*sin(alpha_2);
    A(3,3) = Cperp*l2*(r^2*cos(alpha_1 + beta_1)^2 + l1^2*cos(beta_1)^2 + l2^2/3 + l2*r*cos(alpha_1 + beta_1) + l1*l2*cos(beta_1) + 2*l1*r*cos(alpha_1 + beta_1)*cos(beta_1)) - mu*r^3*torqueCoeff + Cperp*l2*(r^2*cos(alpha_2 + beta_2)^2 + l1^2*cos(beta_2)^2 + l2^2/3 + l2*r*cos(alpha_2 + beta_2) + l1*l2*cos(beta_2) + 2*l1*r*cos(alpha_2 + beta_2)*cos(beta_2)) + (Cperp*l1*(3*r^2*cos(alpha_1)^2 + l1^2 + 3*l1*r*cos(alpha_1)))/3 + (Cperp*l1*(3*r^2*cos(alpha_2)^2 + l1^2 + 3*l1*r*cos(alpha_2)))/3 + Cpar*l1*r^2*sin(alpha_1)^2 + Cpar*l1*r^2*sin(alpha_2)^2 + Cpar*l2*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1))*(r*sin(alpha_1 + beta_1) + l2*sin(beta_1)) + Cpar*l2*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2))*(r*sin(alpha_2 + beta_2) + l2*sin(beta_2));

end