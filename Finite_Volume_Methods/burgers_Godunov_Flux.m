function F = burgers_Godunov_Flux(U_l,U_r)
% returns value of flux function for Godunov's method

if U_l >= U_r % shock
    if (U_l + U_r)/2 > 0
        u_star = U_l;
    else
        u_star = U_r;
    end
else % rarefaction
    if U_l > 0
        u_star = U_l;
    elseif U_r < 0
        u_star = U_r;
    else
        u_star = 0;
    end
end
F = u_star^2/2;