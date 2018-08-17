%% Solution to Burger's equation using Godunov's method
% u_t + u*u_x = 0 on -2 <= x <= 3 and 0 <= t <= 2
% u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0
% from pg. 173 of Noble

%% Edits by KO on 8/17/2018
%       Modified the initial conditions to allow for a smooth sine bump.


clear all; close all; clc
L = 5; % length of x-interval
tf = 1; % length of t-interval
n = 100; % number of cells in x
m = 100; % number of time steps
h = L/n; % cell lengths in x
k = tf/m; % time step size
Q = zeros(n,m); % space-time matrix of Q = cell volumes of function u

%% Building the initial conditions

l_sec = 1:(2*n/5); 
m_sec = (2*n/5 + 1):(3*n/5 - 1);
r_sec = (3*n/5):n;

Q(l_sec,1) = 0; % left initial condition
Q(m_sec,1) = max(0,sin(pi*m_sec/length(m_sec))); % middle initial condition
Q(r_sec,1) =0 ; % right initial condition


a = k/h; % FVM constant
for i = 2:m
    for j = 2:(n-1)
        Q(j,i) = Q(j,i-1) - a*(burgers_Godunov_Flux(Q(j,i-1),Q(j+1,i-1))-...
            burgers_Godunov_Flux(Q(j-1,i-1),Q(j,i-1)));
    end
    plot(Q(:,i))
    pause(.1)
end
figure(1); surf(Q((n/4):(3*n/4),1:m))
title('u[x,0] = 0 if x <= 0, 1 if 0 < x < 1, and 0 if x >= 0 (cell averages)'); 
ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight