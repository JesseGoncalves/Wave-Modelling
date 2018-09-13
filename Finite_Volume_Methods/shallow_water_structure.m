function S = shallow_water_structure(h_l,h_r,u_l,u_r)
% predicts structure of solution to shallow water Riemann problem

g = 9.8; % gravitational constant

% find correct middle states H_m, U_m
F = @(h_m) (h_m < h_l)*(u_l + 2*(sqrt(g*h_l) - sqrt(g*h_m))) + (h_m > h_l)*...
    (u_l - (h_m - h_l)*sqrt(g/2*(1/h_m + 1/h_l))) - ((h_m > h_r)*...
    (u_r - 2*(sqrt(g*h_r) - sqrt(g*h_m))) + (h_m < h_r)*...
    (u_r + (h_m - h_r)*sqrt(g/2*(1/h_m + 1/h_r))));
H_m = fsolve(F,1);
U_m = (H_m < h_l)*(u_l + 2*(sqrt(g*h_l) - sqrt(g*H_m))) + (H_m > h_l)*...
    (u_l - (H_m - h_l)*sqrt(g/2*(1/H_m + 1/h_l)));

% shock speeds (if entropy-satisfying)
s_l = (H_m*U_m - h_l*u_l)/(H_m - h_l);
s_r = (h_r*u_r - H_m*U_m)/(h_r - H_m);

% eigenvalues
lambda1_l = u_l - sqrt(g*h_l);
lambda2_l = u_l + sqrt(g*h_l);
lambda1_m = U_m - sqrt(g*H_m);
lambda2_m = U_m + sqrt(g*H_m);
lambda1_r = u_r - sqrt(g*h_r);
lambda2_r = u_r + sqrt(g*h_r);

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

S = strcat(left,', ',right);
end



