%% Chapter 13
%% Exercise 7
%% (e)
% p(v) = -exp(v)
clear all; close all; clc;
Vl = 1; Ul = 1; Vr = 3; Ur = 4;
V =  linspace(-3,5,100);
figure(1); plot(V,loci7(V,Vl,Ul,1)); hold on;
plot(V,loci7(V,Vr,Ur,2)) % plot loci tangent to r1(1,1),r2(3,4)
F = @(v) [loci7(v,Vl,Ul,1) - loci7(v,Vr,Ur,2)];
Vm1 = fsolve(F,2) % find alternate intermediate state Vm
Um1 = loci7(Vm1,Vl,Ul,1) % find alternate intermediate state Um
figure(2); plot(V,loci7(V,Vr,Ur,1)); hold on;
plot(V,loci7(V,Vl,Ul,2)) % plot loci tangent to r2(1,1),r1(3,4)
F = @(v) [loci7(v,Vl,Ul,2) - loci7(v,Vr,Ur,1)];
Vm2 = fsolve(F,2) % find correct intermediate state Vm
Um2 = loci7(Vm2,Vl,Ul,2) % find correct intermediate state Um

% options = optimset('TolFun',1e-8,'TolX',1e-8,'Display','iter');%'Jacobian','on')
% Vm2 = fsolve(@(X) F,2,options) % find alternate intermediate state Vm

%% (g)
% p(v) = -exp(v)
clear all
close all
clc

vl = 1; ul = 1;
vr = 3; ur = 4;
vm = linspace(-3,5,100);

um1 = loci7(vm,vl,ul,2); % 1-shock connecting ql to qm
um2 = loci7(vm,vr,ur,1); % 2-shock connecting qr to qm

subplot(2,1,2)

sr = (ur - um2)./(vm - vr);

plot(vm,sr)

hold on

plot(vm,-sqrt(exp(vr*ones(size(vm)))),'k')

plot(vm,-sqrt(exp(vm)),'r')

plot(vm,sqrt(exp(vr*ones(size(vm)))),'k--')

plot(vm,sqrt(exp(vm)),'r--')

title('2-shock connecting qr to qm')

legend('speed','r-','m-','r+','m+')

grid on

subplot(2,1,1)

sl = (um1 - ul)./(vl - vm);

plot(vm,sl)

hold on

plot(vm,-sqrt(exp(vl*ones(size(vm)))),'k')

plot(vm,-sqrt(exp(vm)),'r')

plot(vm,sqrt(exp(vl*ones(size(vm)))),'k--')

plot(vm,sqrt(exp(vm)),'r--')

title('1-shock connecting ql to qm')

legend('speed','l-','m-','l+','m+')

%% (e2)
% p(v) = -2v
clear all; close all; clc;
Vl = 1; Ul = 1; Vr = 3; Ur = 4;
V =  linspace(-3,5,100);
figure(1); plot(V,loci7_2(V,Vl,Ul,1)); hold on; grid on
plot(V,loci7_2(V,Vr,Ur,2)) % plot loci tangent to r1(1,1),r2(3,4)
F = @(v) [loci7_2(v,Vl,Ul,1) - loci7_2(v,Vr,Ur,2)];
Vm1 = fsolve(F,2) % find alternate intermediate state Vm
Um1 = loci7_2(Vm1,Vl,Ul,1) % find alternate intermediate state Um
figure(2); plot(V,loci7_2(V,Vr,Ur,1)); hold on; grid on
plot(V,loci7_2(V,Vl,Ul,2)) % plot loci tangent to r2(1,1),r1(3,4)
F = @(v) [loci7_2(v,Vl,Ul,2) - loci7(v,Vr,Ur,1)];
Vm2 = fsolve(F,2) % find correct intermediate state Vm
Um2 = loci7_2(Vm2,Vl,Ul,2) % find correct intermediate state Um

%options = optimset('TolFun',1e-8,'TolX',1e-8,'Display','iter','Jacobian','on');
%Vm2 = fsolve(@(X) F,2,options) % find alternate intermediate state Vm

%% (g2)
% p(v) = -2v
clear all
close all
clc

vl = 1; ul = 1;
vr = 3; ur = 4;
vm = linspace(-3,5,100);

um1 = loci7(vm,vl,ul,1); % 1-shock connecting ql to qm
um2 = loci7(vm,vr,ur,2); % 2-shock connecting qr to qm

subplot(2,1,1)

sr = (ur - um2)./(vm - vr);

plot(vm,sr)

hold on

plot(vm,-sqrt(exp(vr*ones(size(vm)))),'k')

plot(vm,-sqrt(exp(vm)),'r')

plot(vm,sqrt(exp(vr*ones(size(vm)))),'k--')

plot(vm,sqrt(exp(vm)),'r--')

title('2-shock connecting qr to qm')

legend('speed','r-','m-','r+','m+')

grid on

subplot(2,1,2)

sl = (um1 - ul)./(vl - vm);

plot(vm,sl)

hold on

plot(vm,-sqrt(exp(vl*ones(size(vm)))),'k')

plot(vm,-sqrt(exp(vm)),'r')

plot(vm,sqrt(exp(vl*ones(size(vm)))),'k--')

plot(vm,sqrt(exp(vm)),'r--')

title('1-shock connecting ql to qm')

legend('speed','l-','m-','l+','m+')