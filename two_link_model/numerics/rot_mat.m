function R = rot_mat(Phi)
%% Compute the rotation matrix that maps quantities expressed in the swimmer
%% frame to the lab frame.

    R = [cos(Phi), -sin(Phi), 0;
        sin(Phi), cos(Phi), 0;
        0, 0, 1];

end