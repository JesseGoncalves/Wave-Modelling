%% Tests Mesh Size vs. Conservation for Numerical Solution to Burger's equation Riemann Problem
% u_t + u*u_x = 0 on -2 <= x <= 3 and 0 <= t <= 2
% u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0
% from pg. 173 of Noble
clear all; close all; clc
L = 5; % length of x-interval
tf = 2; % length of t-interval
xx = [1 10 50 100]; % vector of sizing factors
%% Time-step Size vs. Conserved Quantity
U_conserve = zeros(length(xx),100*xx(end)-1); % vector to track conservation
for j = 1:length(xx)
    n = 500; % number of x-grid points
    m = 100*xx(j); % vary number of time steps
    h = L/n; % mesh spacing
    k = tf/m; % time step size
    d = k/(2*h); % finite-difference constant
    U = zeros(n,m); % space-time matrix of u
    U(1:(2*n/5),1) = 0; % left initial condition
    U((2*n/5 + 1):(3*n/5 - 1)) = 1; % middle initial condition
    U((3*n/5):n,1) = 0; % right initial condition
    a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
    D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
    for i = 2:m
        E = diag(U(:,i-1));
        U(:,i) = U(:,i-1) - d*E*D*U(:,i-1); % evolve u matrix
        U(:,i) = U(:,i-1) - d/2*E*D*(U(:,i) + U(:,i-1)); % evolve u using improved euler
        U_conserve(j,i-1) = ((sum(U(:,i))-sum(U(:,i-1)))*h/k + 1/2*(U(1,i-1)^2 - U(end,i-1)^2))/m;
    end
end
% matrix of time step sizes and average error at each time-step
U_conserve_timestep = [tf./(100*xx') sum(U_conserve,2)];

%% Spacial grid Size vs. Conserved Quantity
U_conserve = zeros(length(xx),1000); % vector to track conservation
for j = 1:length(xx)
    n = 50*xx(j); % vary number of x-grid points
    m = 1000; % number of time steps
    h = L/n; % mesh spacing
    k = tf/m; % time step size
    d = k/(2*h); % finite-difference constant
    U = zeros(n,m); % space-time matrix of u
    U(1:(2*n/5),1) = 0; % left initial condition
    U((2*n/5 + 1):(3*n/5 - 1)) = 1; % middle initial condition
    U((3*n/5):n,1) = 0; % right initial condition
    a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
    D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
    for i = 2:m
        E = diag(U(:,i-1));
        U(:,i) = U(:,i-1) - d*E*D*U(:,i-1); % evolve u matrix
        U(:,i) = U(:,i-1) - d/2*E*D*(U(:,i) + U(:,i-1)); % evolve u using improved euler
        U_conserve(j,i-1) = ((sum(U(:,i))-sum(U(:,i-1)))*h/k + 1/2*(U(1,i-1)^2 - U(end,i-1)^2))/n;
    end
end
% matrix of grid spacing sizes and average error at each time-step
U_conserve_gridsize = [L./(50*xx') sum(U_conserve,2)];