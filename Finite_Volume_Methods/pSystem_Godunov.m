%% FVM 13.7 Numerical Solution - Godunov Method
% v_t - u_x = 0
% u_t - e^(v)*v_x = 0
% q(x,t) = [v; u]
% solved on -2 <= x <= 2 and 0 <= t <= 1
clear all; close all; clc
L = 4; % length of x-interval
tf = 1; % length of t-interval
n = 400; % number of x-grid points
m = 1000; % number of time steps
h = L/n; % mesh spacing
k = tf/m; % time step size
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v

%% Set ICs
v_left = 1;
v_right = 3;
u_left = 1;
u_right = 4;
V(1:(n/2),1) = v_left; U(1:(n/2),1) = u_left; % left initial conditions
V((n/2+1):n,1) = v_right; U((n/2+1):n,1) = u_right; % right initial conditions
V(1,1:m) = v_left; V(n,1:m) = v_right; U(1,1:m) = u_left; U(n,1:m) = u_right; % keep BCs constant

%% Predict structure
disp(pSystem_structure(v_left,v_right,u_left,u_right))

%% Apply Godunov's method
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

%% Plot solution
% Setting up the x-t grid for the plots
[T,X] = meshgrid(linspace(0,tf,m),linspace(-2,2,n));
t_inds = 1:m;   x_inds = n/4:3*n/4; % setting the indices for viewing.

figure(1); surf(T(x_inds,t_inds),X(x_inds,t_inds),V(x_inds,t_inds)); % show results
title(['v: v- = ' num2str(v_left) ', v+ = ' num2str(v_right)]);ylabel('space');...
    xlabel('time'); shading interp; view(90,-90); axis tight
figure(2); surf(T(x_inds,t_inds),X(x_inds,t_inds),U(x_inds,t_inds));
title(['u: u- = ' num2str(u_left) ', u+ = ' num2str(u_right)]); ylabel('space');...
    xlabel('time'); shading interp; view(90,-90); axis tight

%% Plot loci
figure(3); pSys_plot_loci(V(1,1),V(n,1),U(1,1),U(n,1));

%% Plot conserved quantities
V_conserve = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));
figure(4); plot(V_conserve); hold on; plot(U_conserve); % plot conserved quantities
title('Conserved Quantities'); legend('v','u'); axis([1 m/4 -10 10]);

