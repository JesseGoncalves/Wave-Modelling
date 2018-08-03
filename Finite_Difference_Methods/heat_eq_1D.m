% Solve heat equation on [0,2]x[0,5] with u(0,t)=u(2,t)=0 and u(x,0)=x for
% 0 <= x <= 1 and u(x,0)=2-x for 1 < x <= 2 using Crank-Nicolson
clear all; close all; clc
L = 2; % length of x-interval
tf = 2; % length of t-interval
n = 200; % number of x-grid points
m = 200; % number of time steps
h = L/n; % mesh spacing
k = tf/m; % time step size
r = k/(2*h^2); % Crank-Nicolson constant
U = zeros(n,m); % matrix of heat in space-time
U(1:(n/2),1) = linspace(0,1,n/2)'; % impose initial condition
U((n/2+1):n) = linspace(1,0,n/2)';
M1 = FDcrankNicolson(n,r,0); % create lhs Crank-Nicolson matrix
M2 = FDcrankNicolson(n,r,1); % create rhs Crank-Nicolson matrix
for i = 2:m
    U(2:(n-1),i) = linsolve(M1,M2*U(2:(n-1),i-1));
end
figure(1); mesh(U)
% %% Calculate explicit solution
% syms x t p
% c(p) = int(x*sin(pi*p*x/2),x,[0 1]) + int((2-x)*sin(pi*p*x/2),x,[1 2]);
% Uex(x,t) = symsum(c(p)*sin(pi*p*x/2)*exp(-(pi*p/2)^2*t),p,1,10);
% UexMat = zeros(n,m);
% UexMat(1:(n/2),1) = linspace(0,1,n/2)';
% UexMat((n/2+1):n) = linspace(1,0,n/2)';
% for j = 2:m
%      for  i = 2:(n-1)
%         UexMat(i,j) = double(Uex(2*i/n,5*j/m)); % improve this!!!
%      end
% end
% figure(2); mesh(UexMat)
% %% Compare solutions
% sum(sqrt(sum((U - UexMat).^2,1)))
%% 
x = linspace(0,2,n);
f = x.*(x <= 1) + (2-x).*(x > 1);
c = zeros(1,10);
[X,T] = meshgrid(x,linspace(0,2,m));
Uex = zeros(size(X));
for p = 1:30
    c(p) = trapz(x,f.*sin(p*pi.*x/2));
    Uex = Uex + c(p)*sin(p*pi*X/2).*exp(-(p*pi/2)^2.*T);
end
figure(2); mesh(Uex')
sum(sqrt(sum((U - Uex').^2,1)))

figure(3)
mesh(U-Uex')