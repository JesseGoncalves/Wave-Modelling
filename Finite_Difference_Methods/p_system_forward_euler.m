%% Numerical Solution to p-system Riemann problem
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
d = k/(2*h); % finite-difference constant
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
V(1:(n/2),1) = 1; U(1:(n/2),1) = 1; % left initial conditions
V((n/2+1):n,1) = 3; U((n/2+1):n,1) = 4; % right initial conditions
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v matrix
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u matrix
end
figure(1); surf(V((4*n/9):(5*n/9),1:10)) % show results from first few time steps
title('v: v- = 1, v+ = 3'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(2); surf(U((4*n/9):(5*n/9),1:10))
title('u: u- = 1, u+ = 4'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve1 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve1 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));

%% ICs: q_0(x) = [3; 4] for x<0, [1; 1] for x>=0 
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
V(1:(n/2),1) = 3; U(1:(n/2),1) = 4; % left initial conditions
V((n/2+1):n,1) = 1; U((n/2+1):n,1) = 1; % right initial conditions
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v matrix
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u matrix
end
figure(3); surf(V((4*n/9):(5*n/9),1:10)) % show results from first few time steps
title('v: v- = 3, v+ = 1'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(4); surf(U((4*n/9):(5*n/9),1:10))
title('u: u- = 4, u+ = 1'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve2 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve2 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));

%% ICs: q_0(x) = [1; 1] for x<0, [1.1; 0.9] for x>=0 
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
V(1:(n/2),1) = 1; U(1:(n/2),1) = 1; % left initial conditions
V((n/2+1):n,1) = 1.1; U((n/2+1):n,1) = 0.9; % right initial conditions
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v matrix
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u matrix
end
figure(5); surf(V((4*n/9):(5*n/9),1:10)) % show results from first few time steps
title('v: v- = 1, v+ = 1.1'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(6); surf(U((4*n/9):(5*n/9),1:10))
title('u: u- = 1, u+ = 0.9'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve3 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve3 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));

%% ICS: tanh approximation of q_0(x) = [1; 1] for x<0, [3; 4] for x>=0
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
x = linspace(-2,2,n)';
V(:,1) = tanh(x/.01) + 2; U(:,1) = 1.5*tanh(x/.01) + 2.5; % ICs
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v matrix
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u matrix
end
figure(1); surf(V((4*n/9):(5*n/9),1:10)) % show results from first few time steps
title('v: v- = 1, v+ = 3'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(2); surf(U((4*n/9):(5*n/9),1:10))
title('u: u- = 1, u+ = 4'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve4 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve4 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));
