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
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v using forward euler
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u using forward euler
    V(:,i) = V(:,i-1) + d/2*D*(U(:,i) + U(:,i-1)); % evolve v using improved euler
    U(:,i) = U(:,i-1) + d/2*E*D*(V(:,i) + V(:,i-1)); % evolve u using improved euler
end
figure(1); surf(V((4*n/9):(5*n/9),1:(m/4))) % show results from first few time steps
title('v: v- = 1, v+ = 3'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(2); surf(U((4*n/9):(5*n/9),1:(m/4)))
title('u: u- = 1, u+ = 4'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve1 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve1 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));
figure(3); plot(V_conserve1); hold on; plot(U_conserve1); % plot conserved quantities
title('Conserved Quantities'); legend('v conserve','u conserve'); axis([1 m/4 -10 10]);

%% ICs: q_0(x) = [3; 4] for x<0, [1; 1] for x>=0 
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
V(1:(n/2),1) = 3; U(1:(n/2),1) = 4; % left initial conditions
V((n/2+1):n,1) = 1; U((n/2+1):n,1) = 1; % right initial conditions
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v using forward euler
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u using forward euler
    V(:,i) = V(:,i-1) + d/2*D*(U(:,i) + U(:,i-1)); % evolve v using improved euler
    U(:,i) = U(:,i-1) + d/2*E*D*(V(:,i) + V(:,i-1)); % evolve u using improved euler
end
figure(4); surf(V((4*n/9):(5*n/9),1:(m/4))) % show results from first few time steps
title('v: v- = 3, v+ = 1'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(5); surf(U((4*n/9):(5*n/9),1:(m/4)))
title('u: u- = 4, u+ = 1'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve2 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve2 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));
figure(6); plot(V_conserve2); hold on; plot(U_conserve2); % plot conserved quantities
title('Conserved Quantities'); legend('v conserve','u conserve'); axis([1 m/4 -10 10]);

%% ICs: q_0(x) = [1; 1] for x<0, [1.1; 0.9] for x>=0 
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
V(1:(n/2),1) = 1; U(1:(n/2),1) = 1; % left initial conditions
V((n/2+1):n,1) = 1.1; U((n/2+1):n,1) = 0.9; % right initial conditions
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v using forward euler
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u using forward euler
    V(:,i) = V(:,i-1) + d/2*D*(U(:,i) + U(:,i-1)); % evolve v using improved euler
    U(:,i) = U(:,i-1) + d/2*E*D*(V(:,i) + V(:,i-1)); % evolve u using improved euler
end
figure(7); surf(V((4*n/9):(5*n/9),1:(m/4))) % show results from first few time steps
title('v: v- = 1, v+ = 1.1'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(8); surf(U((4*n/9):(5*n/9),1:(m/4)))
title('u: u- = 1, u+ = 0.9'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve3 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve3 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));
figure(9); plot(V_conserve3); hold on; plot(U_conserve3); % plot conserved quantities
title('Conserved Quantities'); legend('v conserve','u conserve'); axis([1 m/4 -10 10]);

%% ICS: tanh approximation of q_0(x) = [1; 1] for x<0, [3; 4] for x>=0
U = zeros(n,m); % space-time matrix of u
V = zeros(n,m); % space-time matrix of v
x = linspace(-2,2,n)';
V(:,1) = tanh(x/.01) + 2; U(:,1) = 1.5*tanh(x/.01) + 2.5; % ICs
a = ones(n-1,1); D = diag(a,1) + diag(-a,-1); % finite difference matrix
D(1,:) = zeros(1,size(D,1)); D(end,:) = zeros(1,size(D,1)); % constant BCs
for i = 2:m
    V(:,i) = V(:,i-1) + d*D*U(:,i-1); % evolve v using forward euler
    E = expm(diag(V(:,i-1)));
    U(:,i) = U(:,i-1) + d*E*D*V(:,i-1); % evolve u using forward euler
    V(:,i) = V(:,i-1) + d/2*D*(U(:,i) + U(:,i-1)); % evolve v using improved euler
    U(:,i) = U(:,i-1) + d/2*E*D*(V(:,i) + V(:,i-1)); % evolve u using improved euler
end
figure(10); surf(V((4*n/9):(5*n/9),1:(m/4))) % show results from first few time steps
title('v: tanh approx v- = 1, v+ = 3'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
figure(11); surf(U((4*n/9):(5*n/9),1:(m/4)))
title('u: tanh approx u- = 1, u+ = 4'); ylabel('space'); xlabel('time'); shading interp; view(90,-90); axis tight
V_conserve4 = diff(sum(V))*h/k + U(1,1) - U (end,1); % track conservation: should be zero
U_conserve4 = diff(sum(U))*h/k + exp(V(1,1)) - exp(V(end,1));
figure(12); plot(V_conserve1); hold on; plot(U_conserve1); % plot conserved quantities
title('Conserved Quantities'); legend('v conserve','u conserve'); axis([1 m/4 -10 10]);
