function S = pSystem_structure(v_l,v_r,u_l,u_r)
% predicts structure of solution to p-system Riemann problem [p(v)=-exp(v)]

% find correct middle states H_m, U_m
F = @(v_m) (v_m > v_l)*(u_l + 2*(exp(v_m/2) - exp(v_l/2))) + (v_m < v_l)*...
    (u_l + sqrt((v_m - v_l)*(exp(v_m) - exp(v_l)))) - (v_m > v_r)*...
    (u_r - 2*(exp(v_r/2) - exp(v_m/2))) + (v_m < v_r)*...
    (u_r - sqrt((v_r - v_m)*(exp(v_r) - exp(v_m))));
V_m = fsolve(F,1);
U_m = (V_m > v_l)*(u_l + 2*(exp(V_m/2) - exp(v_l/2))) + (V_m < v_l)*...
    (u_l + sqrt((V_m - v_l)*(exp(V_m) - exp(v_l))));

% shock speeds (if entropy-satisfying)
s_l = (U_m - u_l)/(v_l - V_m);
s_r = (u_r - U_m)/(V_m - v_r);

% eigenvalues
lambda1_l = -exp(v_l);
lambda2_l = exp(v_l);
lambda1_m = -exp(V_m);
lambda2_m = exp(V_m);
lambda1_r = -exp(v_r);
lambda2_r = exp(v_r);

% check entropy conditions
if lambda1_l > lambda1_m
    left = '1-shock';
else
    left = '1-rarefaction';
end

if lambda2_r < lambda2_m
    right = '2-shock';
else
    right = '2-rarefaction';
end

S = strcat(left,', ',right,' V_m = ',num2str(V_m));
end



