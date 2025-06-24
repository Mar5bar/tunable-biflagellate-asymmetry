function Fn = compute_contact_forces(or,position,shape,velocities,shape_velocities,params)

alpha_1 = shape(1);
alpha_2 = shape(2);
beta_1 = shape(3);
beta_2 = shape(4);

phi1 = position(1);
phi2 = position(2);

% Velocity in the swimmer frame.
R = rot_mat(-or);
U = R(1:2,1:2)*[velocities(1);velocities(2)];
ux = U(1);
uy = U(2);
wz = velocities(3);

alpha_1dot = shape_velocities(1);
alpha_2dot = shape_velocities(2);
beta_1dot = shape_velocities(3);
beta_2dot = shape_velocities(4);

l1 = params.l1;
l2 = params.l2;
Cperp = params.Cperp;
Cpar = params.Cpar;
r = params.r;

% x and y components of first filament force
Fn(1) = (Cperp*l2*sin(alpha_1 + beta_1 - phi1)*(alpha_1dot*l2 + beta_1dot*l2 + 2*alpha_1dot*l1*cos(beta_1)))/2 - l2*(Cpar*cos(alpha_1 + beta_1 - phi1)*(ux*cos(alpha_1 + beta_1 - phi1) - uy*sin(alpha_1 + beta_1 - phi1)) + Cperp*sin(alpha_1 + beta_1 - phi1)*(uy*cos(alpha_1 + beta_1 - phi1) + ux*sin(alpha_1 + beta_1 - phi1))) - l2*wz*(Cperp*sin(alpha_1 + beta_1 - phi1)*(l2/2 + r*cos(alpha_1 + beta_1) + l1*cos(beta_1)) - Cpar*cos(alpha_1 + beta_1 - phi1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1))) - l1*wz*((Cperp*l1*sin(alpha_1 - phi1))/2 - Cpar*r*sin(alpha_1)*cos(alpha_1 - phi1) + Cperp*r*sin(alpha_1 - phi1)*cos(alpha_1)) - l1*(Cpar*ux*cos(alpha_1 - phi1)^2 - (Cpar*uy*sin(2*alpha_1 - 2*phi1))/2 + (Cperp*uy*sin(2*alpha_1 - 2*phi1))/2 + Cperp*ux*sin(alpha_1 - phi1)^2) + (Cperp*alpha_1dot*l1^2*sin(alpha_1 - phi1))/2 + Cpar*alpha_1dot*l1*l2*(sin(alpha_1 - phi1)/2 - sin(alpha_1 + 2*beta_1 - phi1)/2);
Fn(2) = (Cperp*l2*cos(alpha_1 + beta_1 - phi1)*(alpha_1dot*l2 + beta_1dot*l2 + 2*alpha_1dot*l1*cos(beta_1)))/2 - l2*(Cperp*cos(alpha_1 + beta_1 - phi1)*(uy*cos(alpha_1 + beta_1 - phi1) + ux*sin(alpha_1 + beta_1 - phi1)) - Cpar*sin(alpha_1 + beta_1 - phi1)*(ux*cos(alpha_1 + beta_1 - phi1) - uy*sin(alpha_1 + beta_1 - phi1))) - l2*wz*(Cperp*cos(alpha_1 + beta_1 - phi1)*(l2/2 + r*cos(alpha_1 + beta_1) + l1*cos(beta_1)) + Cpar*sin(alpha_1 + beta_1 - phi1)*(r*sin(alpha_1 + beta_1) + l1*sin(beta_1))) - l1*wz*((Cperp*l1*cos(alpha_1 - phi1))/2 + Cperp*r*cos(alpha_1)*cos(alpha_1 - phi1) + Cpar*r*sin(alpha_1 - phi1)*sin(alpha_1)) - l1*(Cperp*uy*cos(alpha_1 - phi1)^2 - (Cpar*ux*sin(2*alpha_1 - 2*phi1))/2 + (Cperp*ux*sin(2*alpha_1 - 2*phi1))/2 + Cpar*uy*sin(alpha_1 - phi1)^2) + (Cperp*alpha_1dot*l1^2*cos(alpha_1 - phi1))/2 - Cpar*alpha_1dot*l1*l2*(cos(alpha_1 + 2*beta_1 - phi1)/2 - cos(alpha_1 - phi1)/2);
% x and y components of second filament force
Fn(3) = (Cperp*l2*sin(alpha_2 + beta_2 - phi2)*(alpha_2dot*l2 + beta_2dot*l2 + 2*alpha_2dot*l1*cos(beta_2)))/2 - l2*(Cpar*cos(alpha_2 + beta_2 - phi2)*(ux*cos(alpha_2 + beta_2 - phi2) - uy*sin(alpha_2 + beta_2 - phi2)) + Cperp*sin(alpha_2 + beta_2 - phi2)*(uy*cos(alpha_2 + beta_2 - phi2) + ux*sin(alpha_2 + beta_2 - phi2))) - l2*wz*(Cperp*sin(alpha_2 + beta_2 - phi2)*(l2/2 + r*cos(alpha_2 + beta_2) + l1*cos(beta_2)) - Cpar*cos(alpha_2 + beta_2 - phi2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2))) - l1*wz*((Cperp*l1*sin(alpha_2 - phi2))/2 - Cpar*r*sin(alpha_2)*cos(alpha_2 - phi2) + Cperp*r*sin(alpha_2 - phi2)*cos(alpha_2)) - l1*(Cpar*ux*cos(alpha_2 - phi2)^2 - (Cpar*uy*sin(2*alpha_2 - 2*phi2))/2 + (Cperp*uy*sin(2*alpha_2 - 2*phi2))/2 + Cperp*ux*sin(alpha_2 - phi2)^2) + (Cperp*alpha_2dot*l1^2*sin(alpha_2 - phi2))/2 + Cpar*alpha_2dot*l1*l2*(sin(alpha_2 - phi2)/2 - sin(alpha_2 + 2*beta_2 - phi2)/2);
Fn(4) = (Cperp*l2*cos(alpha_2 + beta_2 - phi2)*(alpha_2dot*l2 + beta_2dot*l2 + 2*alpha_2dot*l1*cos(beta_2)))/2 - l2*(Cperp*cos(alpha_2 + beta_2 - phi2)*(uy*cos(alpha_2 + beta_2 - phi2) + ux*sin(alpha_2 + beta_2 - phi2)) - Cpar*sin(alpha_2 + beta_2 - phi2)*(ux*cos(alpha_2 + beta_2 - phi2) - uy*sin(alpha_2 + beta_2 - phi2))) - l2*wz*(Cperp*cos(alpha_2 + beta_2 - phi2)*(l2/2 + r*cos(alpha_2 + beta_2) + l1*cos(beta_2)) + Cpar*sin(alpha_2 + beta_2 - phi2)*(r*sin(alpha_2 + beta_2) + l1*sin(beta_2))) - l1*wz*((Cperp*l1*cos(alpha_2 - phi2))/2 + Cperp*r*cos(alpha_2)*cos(alpha_2 - phi2) + Cpar*r*sin(alpha_2 - phi2)*sin(alpha_2)) - l1*(Cperp*uy*cos(alpha_2 - phi2)^2 - (Cpar*ux*sin(2*alpha_2 - 2*phi2))/2 + (Cperp*ux*sin(2*alpha_2 - 2*phi2))/2 + Cpar*uy*sin(alpha_2 - phi2)^2) + (Cperp*alpha_2dot*l1^2*cos(alpha_2 - phi2))/2 - Cpar*alpha_2dot*l1*l2*(cos(alpha_2 + 2*beta_2 - phi2)/2 - cos(alpha_2 - phi2)/2);


end