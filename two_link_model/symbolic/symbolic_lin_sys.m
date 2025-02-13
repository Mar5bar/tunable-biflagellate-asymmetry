addpath(genpath('.'))

% Fixed parameters.
syms nx1 ny1 nz1 tx1 ty1 tz1 bx1 by1 bz1 real
syms nx2 ny2 nz2 tx2 ty2 tz2 bx2 by2 bz2 real
syms l1 l2 real
syms Cpar Cperp real
syms r real
syms phi1 phi2 real
syms forceCoeff torqueCoeff real
syms mu real

nx = [nx1; nx2];
ny = [ny1; ny2];
nz = [nz1; nz2];
tx = [tx1; tx2];
ty = [ty1; ty2];
tz = [tz1; tz2];
bx = [bx1; bx2];
by = [by1; by2];
bz = [bz1; bz2];

phi = [phi1; phi2];

% Controlled angles.
syms alpha_1 alpha_2 beta_1 beta_2 real
syms alpha_1dot alpha_2dot beta_1dot beta_2dot real
alpha = [alpha_1; alpha_2];
beta = [beta_1; beta_2];
alphadot = [alpha_1dot; alpha_2dot];
betadot = [beta_1dot; beta_2dot];

% Unknowns.
syms ux uy wz real
u = [ux; uy; 0];
w = [0; 0; wz];

% Lab Euler angles.
syms Phi

% Rotation matrix.
R = [cos(Phi), -sin(Phi), 0;
    sin(Phi), cos(Phi), 0;
    0, 0, 1];


Fs = [];
Gs = [];
Ms = [];
Ns = [];
for i = 1 : 2
    ni = [nx(i); ny(i); nz(i)];
    ti = [tx(i); ty(i); tz(i)];
    bi = [bx(i); by(i); bz(i)];
    ri = r*ni;
    ppar = cos(alpha(i))*ni + sin(alpha(i))*ti;
    pperp = sin(alpha(i))*ni - cos(alpha(i))*ti;
    qpar = cos(alpha(i)+beta(i))*ni + sin(alpha(i)+beta(i))*ti;
    qperp = sin(alpha(i)+beta(i))*ni - cos(alpha(i)+beta(i))*ti;

    % rcrossnCoeff = simplify(dot(cross(ri,ni), bi));
    % rcrosstCoeff = simplify(dot(cross(ri,ti), bi));
    rcrossnCoeff = simplify(dot(cross(ri,ni), bi) / dot(cross(ti,ni),bi));
    rcrosstCoeff = simplify(dot(cross(ri,ti), bi) / dot(cross(ti,ni),bi));

    rcrosspparCoeff = rcrossnCoeff*cos(alpha(i)) + rcrosstCoeff*sin(alpha(i));
    rcrosspperpCoeff = rcrossnCoeff*sin(alpha(i)) - rcrosstCoeff*cos(alpha(i));
    rcrossqparCoeff = rcrossnCoeff*cos(alpha(i)+beta(i)) + rcrosstCoeff*sin(alpha(i)+beta(i));
    rcrossqperpCoeff = rcrossnCoeff*sin(alpha(i)+beta(i)) - rcrosstCoeff*cos(alpha(i)+beta(i));

    wdotb = dot(w,bi);
    udotppar = dot(u,ppar);
    udotpperp = dot(u,pperp);
    udotqpar = dot(u,qpar);
    udotqperp = dot(u,qperp);

    % Force from first segment.
    Fs = [Fs, l1 * (Cpar*udotppar*ppar + ...
            Cperp*udotpperp*pperp ... 
            ) + ...
            l1 * wdotb * ...
            (Cpar*rcrosspparCoeff*ppar + ...
             Cperp*(rcrosspperpCoeff + l1/2)*pperp) - ...
            Cperp*l1^2/2*alphadot(i)*pperp];

    % Force from second segment.
    Gs = [Gs, l2 * (Cpar*udotqpar*qpar + ...
                Cperp*udotqperp*qperp ...
                ) + ...
                l2 * wdotb * ...
                (Cpar*(rcrossqparCoeff - l1*sin(beta(i)))*qpar + ...
                    Cperp * (rcrossqperpCoeff + l1*cos(beta(i)) + l2/2) * qperp) - ...
                Cperp*(l1*l2*alphadot(i)*cos(beta(i)) + l2^2/2*(alphadot(i) + betadot(i)))*qperp + ...
                Cpar*l1*l2*alphadot(i)*sin(beta(i))*qpar];

    % Torque from first segment.
    Ms = [Ms, bi*(...
        Cperp*l1*((udotpperp + rcrosspperpCoeff*wdotb)*(rcrosspperpCoeff + l1/2) + (l1/2*rcrosspperpCoeff + l1^2/3)*(wdotb - alphadot(i))) + ...
        +rcrosspparCoeff*Cpar*l1*(udotppar + wdotb*rcrosspparCoeff))];

    % Torque from second segment.
    Ns = [Ns, bi*(...
        (Cpar*l1*l2*alphadot(i)*sin(beta(i))*(-l1*sin(beta(i))+rcrossqparCoeff) +...
         Cperp*(-l2^2/2*l1*alphadot(i)*cos(beta(i)) - l2^3/3*(alphadot(i)+betadot(i)) - ...
            (l1*cos(beta(i))+rcrossqperpCoeff)*(l1*l2*alphadot(i)*cos(beta(i))+l2^2/2*(alphadot(i)+betadot(i))))) + ...
        Cperp*l2*(((l1*cos(beta(i)) + rcrossqperpCoeff)*wdotb + udotqperp)*(l2/2+l1*cos(beta(i))+rcrossqperpCoeff) + ...
            wdotb*(l2/2*rcrossqperpCoeff + l1*l2/2*cos(beta(i)) + l2^2/3)) + ...
        l2*Cpar*(l1*sin(beta(i))-rcrossqparCoeff)*(udotqpar + wdotb*(-rcrossqparCoeff + l2*sin(beta(i)))))];
end

% Force and torque on the oblate ellipsoidal body.
bodyForce = -mu * forceCoeff * [1;1;1] * r .* u;
bodyTorque = -mu * torqueCoeff * [1;1;1] * r^3 .* w;

forceBalance = sum(Fs + Gs, 2) + bodyForce;
torqueBalance = sum(Ms + Ns, 2) + bodyTorque;

linSys = sym(zeros(3,3));
RHSMat = sym(zeros(3,4));
% Assign the force balance to linSys and RHSMat.
for j = 1 : 2
    temp = coeffs(forceBalance(j),ux); linSys(j,1) = temp(end);
    temp = coeffs(forceBalance(j),uy); linSys(j,2) = temp(end);
    temp = coeffs(forceBalance(j),wz); linSys(j,3) = temp(end);

    for i = 1 : 2
        temp = coeffs(forceBalance(j),alphadot(i)); RHSMat(j,i) = temp(end);
        temp = coeffs(forceBalance(j),betadot(i)); RHSMat(j,i+2) = temp(end);
    end
end
% Assign the torque balance to linSys and RHSMat.
j = 3;
temp = coeffs(torqueBalance(j),ux); linSys(j,1) = temp(end);
temp = coeffs(torqueBalance(j),uy); linSys(j,2) = temp(end);
temp = coeffs(torqueBalance(j),wz); linSys(j,3) = temp(end);
for i = 1 : 2
    temp = coeffs(torqueBalance(j),alphadot(i)); RHSMat(j,i) = temp(end);
    temp = coeffs(torqueBalance(j),betadot(i)); RHSMat(j,i+2) = temp(end);
end
RHSMat = - RHSMat;


% Sub in theta and phi.

expr = linSys;
subs_theta_phi;
linSys = simplify(expr);

expr = RHSMat;
subs_theta_phi;
RHSMat = simplify(expr);

disp(linSys)
disp(RHSMat)
save('symbolic_calculation_general.mat')
stop

% Form the control matrix M_swim, where [u;w] = M_swim*[alphadot;betadot], with u,w in the swimmer frame.
% controlMat_swim = simplify(inv(linSys) * RHSMat);

% To find the evolution of lab frame coordinates, we want to decompose 
% linSys = linSys_lab * [R',0;0,R'].
% Note that linSys_lab = linSys * [R,0;0,R].
% linSys_lab = linSys * blkdiag(R,R);

% We have linSys_lab * [U_lab; Omega_lab] = RHSMat * [alphadot;betadot].
% So, form the control matrix M_lab, where [u;w] = M_lab*[alphadot;betadot], with u,w in the lab frame.
% controlMat_lab = simplify(inv(linSys_lab) * RHSMat);

%% Special cases.
% We'll first specify some geometry: phi_i = i*pi/2, theta_i = theta_0.
% We'll overwrite the quantities above.
sub = @(expr) simplify(subs(expr, [theta1,theta2,theta3,theta4,phi1,phi2,phi3,phi4], [theta0,theta0,theta0,theta0,0,pi/2,pi,3*pi/2]));
linSys = sub(linSys);
linSys_lab = linSys * blkdiag(R,R);
% controlMat_swim = simplify(inv(subs(linSys, [theta1,theta2,theta3,theta4,phi1,phi2,phi3,phi4], [theta0,theta0,theta0,theta0,0,pi/2,pi,3*pi/2])) * subs(RHSMat, [theta1,theta2,theta3,theta4,phi1,phi2,phi3,phi4], [theta0,theta0,theta0,theta0,0,pi/2,pi,3*pi/2]))
% controlMat_lab = simplify(inv(subs(linSys * blkdiag(R,R), [theta1,theta2,theta3,theta4,phi1,phi2,phi3,phi4], [theta0,theta0,theta0,theta0,0,pi/2,pi,3*pi/2])) * subs(RHSMat, [theta1,theta2,theta3,theta4,phi1,phi2,phi3,phi4], [theta0,theta0,theta0,theta0,0,pi/2,pi,3*pi/2]))

%% Pronk.
% alpha_1 = alpha_2 = alpha_3 = alpha_4 == alpha0, beta_1 = beta_2 = beta_3 = beta_4 == beta0.
pronk = struct();
sub = @(expr) simplify(subs(expr, [alpha_1,alpha_2,alpha_3,alpha_4,beta_1,beta_2,beta_3,beta_4], [alpha0,alpha0,alpha0,alpha0,beta0,beta0,beta0,beta0]));
pronk.controlMat_swim = simplify(linsolve(sub(linSys),sub(RHSMat)))
pronk.controlMat_lab = simplify(linsolve(sub(linSys_lab),sub(RHSMat)))

disp(pronk)

%% Trot.
% alpha_1 = alpha_3, alpha_2 = alpha_4, beta_1 = beta_3, beta_2 = beta_4.
% Note that this also needs some dynamic conditions.
trot = struct();
sub = @(expr) simplify(subs(expr, [alpha_3,alpha_4,beta_3,beta_4], [alpha_1,alpha_2,beta_1,beta_2]));
trot.controlMat_swim = simplify(linsolve(sub(linSys),sub(RHSMat)))
trot.controlMat_lab = simplify(linsolve(sub(linSys_lab),sub(RHSMat)))

disp(trot)