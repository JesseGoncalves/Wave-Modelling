function F_u = pSystem_U_Godunov_Flux(V_l,V_r,U_l,U_r,max_lambda)
% returns value of flux function in Godunov's Method for second eq in
% p-system

F_u = -.5*(exp(V_l) + exp(V_r) + max_lambda*(U_r - U_l));
