function F_v = pSystem_V_Godunov_Flux(V_l,V_r,U_l,U_r,max_lambda)
% returns value of flux function in Godunov's Method for first eq in
% p-system

F_v = -.5*(U_l + U_r + max_lambda*(V_r - V_l));