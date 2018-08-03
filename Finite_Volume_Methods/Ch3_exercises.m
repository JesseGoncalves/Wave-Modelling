%% Finite Volume Methods
%% Exercise 3.2
clear; clc
A = [0 4; 1 0];
q_l = [0; 1];
q_r = [1; 1];
R = zeros(size(A,1));
eigVals = zeros(size(A,1));
[R eigVals] = eig(A);
R_inv = inv(R);
w_l = R_inv*q_l;
w_r = R_inv*q_r;
x = linspace(-3,3,100);
t = 1;

q=(w_l(1)*R(:,1)+w_l(2)*R(:,2)).*meshgrid(x>eigVals(1,1)*t,[1 1])%...
%    +(eigVals(1,1)*t>x).*(x>eigVals(2,2)*t)*(w_l(1).*R(:,1)+w_r(2).*R(:,2))...
%    +(w_r(1).*R(:,1)+w_r(2).*R(:,2))*(x>eigVals(2,2)*t);
plot(x,q)
return

q = piecewise(x>eigVals(1,1)*t, w_l(1).*R(:,1)+w_l(2).*R(:,2),...
    eigVals(1,1)*t>x>eigVals(2,2)*t, w_l(1).*R(:,1)+w_r(2).*R(:,2),...
    x>eigVals(2,2)*t, w_r(1).*R(:,1)+w_r(2).*R(:,2))
% plot?????

