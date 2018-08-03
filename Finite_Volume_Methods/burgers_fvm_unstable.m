%% Finite Volume Method Solution to Burger's equation Riemann Problem--Unstable!
% u_t + u*u_x = 0 on -2 <= x <= 3 and 0 <= t <= 2
% u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0
% from pg. 173 of Noble
clear all; close all; clc
L = 5; % length of x-interval
tf = 2; % length of t-interval
n = 500; % number of cells in x
m = 1000; % number of time steps
h = L/n; % cell lengths in x
k = tf/m; % time step size
Q = zeros(n,m); % space-time matrix of Q = cell volumes of function u
Q(1:(2*n/5),1) = 0; % left initial condition
Q((2*n/5 + 1):(3*n/5 - 1),1) = h; % middle initial condition
Q((3*n/5):n,1) = 0; % right initial condition
d = k/h; % fvm constant
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    Q(:,i) = Q(:,i-1) - d*D*(Q(:,i-1)).^2;
end
figure(1); surf(Q((n/4):(3*n/4),1:1000))
title('u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0 (cell averages)'); 
ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight