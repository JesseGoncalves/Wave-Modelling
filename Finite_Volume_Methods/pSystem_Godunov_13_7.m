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
a = k/h; % FVM constant
for i = 2:m
    for j = 2:(n-1)
        V(i,j) = V(i,j-1) - a*(
    end
end
