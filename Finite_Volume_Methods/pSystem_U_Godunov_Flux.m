function F_u = pSystem_U_Godunov_Flux(V_l,V_r,U_l,U_r)
% returns value of flux function in Godunov's Method for second eq in
% p-system

s = (exp(V_r) - exp(V_l))/(U_l - U_r); % shock speed
if U_l > U_r % shock
    if s > 0
        u_star = U_l;
    else
        u_star = U_r;
    end
else % rarefaction
    if exp(V_l) > 0
        u_star = U_l;
    elseif exp(V_l) < 0 < exp(V_r)
        u_star = U_l - 2*(exp(V_r/2) - exp(V_l/2));
    else
        u_star = U_r;
    end
end
F_u = 
