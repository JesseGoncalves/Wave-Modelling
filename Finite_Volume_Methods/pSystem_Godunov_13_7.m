%% FVM 13.7 Numerical Solution - Godunov Method
% v_t - u_x = 0 on -2 <= x <= 2 and 0 <= t <= 1
% u_t - e^(v)*v_x = 0
% q(x,t) = [v; u]
% q_0(x) = [1; 1] for x<0, [3; 4] for x>=0 
clear all; close all; clc
L = 4; % length of x-interval
tf = 1; % length of t-interval
n = 400; % number of x-grid points
m = 1000; % number of time steps
h = L/n; % mesh spacing
k = tf/m; % time step size
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
V(1:(n/2),1) = 1; U(1:(n/2),1) = 1; % left initial conditions
V((n/2+1):n,1) = 3; U((n/2+1):n,1) = 4; % right initial conditions
V(1,1:m) = 1; V(n,1:m) = 3; U(1,1:m) = 1; U(n,1:m) = 4; % keep BCs constant
a = k/(2*h); % FVM constant
for j = 2:m
    for i = 2:(n-1)
        max_lambda1 = max([exp(V(i,j-1)) exp(V(i+1,j-1))]);
        max_lambda2 = max([exp(V(i-1,j-1)) exp(V(i,j-1))]);
        V(i,j) = V(i,j-1) - a*(pSystem_V_Godunov_Flux(V(i,j-1),...
            V(i+1,j-1),U(i,j-1),U(i+1,j-1),max_lambda1) - ...
            pSystem_V_Godunov_Flux(V(i-1,j-1),V(i,j-1),U(i-1,j-1),...
            U(i,j-1),max_lambda2)); % evolve V
        U(i,j) = U(i,j-1) - a*(pSystem_U_Godunov_Flux(V(i,j-1),...
            V(i+1,j-1),U(i,j-1),U(i+1,j-1),max_lambda1) - ...
            pSystem_U_Godunov_Flux(V(i-1,j-1),V(i,j-1),U(i-1,j-1),...
            U(i,j-1),max_lambda2)); % evolve U
    end
end

figure(1); surf(V(n/4:3*n/4,1:m)); % show results
title('v: v- = 1, v+ = 3'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(2); surf(U(n/4:3*n/4,1:m));
title('u: u- = 1, u+ = 4'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight