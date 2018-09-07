function F_h = shallow_water_HU_Godunov_flux(H_l,H_r,HU_l,HU_r,max_lambda)
% returns value of flux function in Godunov's Method for second eq in
% shallow water system

g = 9.8; % gravitational constant
F_h = HU_l^2/H_l + 0.5*g*H_l^2 + HU_r^2/H_r + 0.5*g*H_r^2 - ...
    max_lambda*(HU_r - HU_l);