%% Godunov's method for the shallow water equations
% h_t + (hu)_x = 0
% (hu)_t + (hu^2 + 1/2*gh^2)_x = 0
% q(x,t) = [h; hu]
% solved on -2 <= x <= 2 and 0 <= t <= 1
clear all; close all; clc
g = 9.8; % gravitational constant
L = 4; % length of x-interval
tf = 1; % length of t-interval
n = 400; % number of x-grid points
m = 2000; % number of time steps
h = L/n; % mesh spacing
k = tf/m; % time step size
H = zeros(n,m); % space-time matrix of h
HU = zeros(n,m); % space-time matrix of hu

%% Set ICs
h_left = 1; % must be > 0
h_right = 1; % must be > 0
u_left = -2;
u_right = 2;
H(1:(n/2),1) = h_left; HU(1:(n/2),1) = h_left*u_left; % left initial conditions
H((n/2+1):n,1) = h_right; HU((n/2+1):n,1) = h_right*u_right; % right initial conditions
H(1,1:m) = h_left; H(n,1:m) = h_right; HU(1,1:m) = h_left*u_left;...
    HU(n,1:m) = h_right*u_right; % keep BCs constant

%% Predict structure
disp(shallow_water_structure(h_left,h_right,u_left,u_right))

%% Apply Godunov's method
a = k/(2*h); % FVM constant
for j = 2:m
    for i = 2:(n-1)
        max_lambda1 = max([HU(i,j-1)/H(i,j-1)+sqrt(g*H(i,j-1))...
            HU(i+1,j-1)/H(i+1,j-1)+sqrt(g*H(i+1,j-1))]);
        max_lambda2 = max([HU(i-1,j-1)/H(i-1,j-1)+sqrt(g*H(i-1,j-1))...
            HU(i,j-1)/H(i,j-1)+sqrt(g*H(i,j-1))]);
        H(i,j) = H(i,j-1) - a*(shallow_water_H_Godunov_flux(H(i,j-1),...
            H(i+1,j-1),HU(i,j-1),HU(i+1,j-1),max_lambda1) - ...
            shallow_water_H_Godunov_flux(H(i-1,j-1),H(i,j-1),HU(i-1,j-1),...
            HU(i,j-1),max_lambda2)); % evolve H
        HU(i,j) = HU(i,j-1) - a*(shallow_water_HU_Godunov_flux(H(i,j-1),...
            H(i+1,j-1),HU(i,j-1),HU(i+1,j-1),max_lambda1) - ...
            shallow_water_HU_Godunov_flux(H(i-1,j-1),H(i,j-1),HU(i-1,j-1),...
            HU(i,j-1),max_lambda2)); % evolve HU
    end
end

%% Plot solution

% Setting up the x-t grid for the plots
[T,X] = meshgrid(linspace(0,tf,m),linspace(-2,2,n));
t_inds = 1:m;   x_inds = n/4:3*n/4; % setting the indices for viewing.

figure(1); surf(T(x_inds,t_inds),X(x_inds,t_inds),H(x_inds,t_inds)); % show results
title(['h: h- = ' num2str(h_left) ', h+ = ' num2str(h_right)]);...
    ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(2); surf(T(x_inds,t_inds),X(x_inds,t_inds),HU(x_inds,t_inds));
title(['hu: hu- = ' num2str(h_left*u_left) ', hu+ = ' num2str(h_right*u_right)]);...
    ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(3); surf(T(x_inds,t_inds),X(x_inds,t_inds),HU(x_inds,t_inds)./H(x_inds,t_inds));
title(['u: u- = ' num2str(u_left) ', u+ = ' num2str(u_right)]); ylabel('space');...
    xlabel('time'); shading interp; view(90,-90); axis tight
