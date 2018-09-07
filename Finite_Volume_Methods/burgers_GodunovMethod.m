%% Solution to Burger's equation using Godunov's method
% u_t + u*u_x = 0 on -2 <= x <= 3 and 0 <= t <= 2
% u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0
% from pg. 173 of Noble
clear all; close all; clc
L = 5; % length of x-interval
tf = 2; % length of t-interval
n = 1000; % number of cells in x
m = 500; % number of time steps
h = L/n; % cell lengths in x
k = tf/m; % time step size
Q = zeros(n,m); % space-time matrix of Q = cell volumes of function u
Q(1:(2*n/5),1) = 0; % left initial condition
Q((2*n/5 + 1):(3*n/5 - 1),1) = 1; % middle initial condition
Q((3*n/5):n,1) = 0; % right initial condition
a = k/h; % FVM constant
for i = 2:m
    for j = 2:(n-1)
        Q(j,i) = Q(j,i-1) - a*(burgers_Godunov_Flux(Q(j,i-1),Q(j+1,i-1))-...
            burgers_Godunov_Flux(Q(j-1,i-1),Q(j,i-1)));
    end
end
figure(1); surf(Q((n/4):(3*n/4),1:m))
title('u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0 (cell averages)'); 
ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight