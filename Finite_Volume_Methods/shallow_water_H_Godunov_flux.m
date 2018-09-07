function F_h = shallow_water_H_Godunov_flux(H_l,H_r,HU_l,HU_r,max_lambda)
% returns value of flux function in Godunov's Method for first eq in
% shallow water system

F_h = HU_l + HU_r - max_lambda*(H_r - H_l);