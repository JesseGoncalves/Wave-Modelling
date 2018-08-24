function F_v = pSystem_V_Godunov_Flux(V_l,V_r,U_l,U_r)
% returns value of flux function in Godunov's Method for first eq in
% p-system

s = (U_r - U_l)/(V_l - V_r); % shock speed
if s > 0
    v_star = V_l;
else
    v_star = V_r;
end
F_v = -v_star;