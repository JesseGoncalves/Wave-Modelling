%% Finite Volume Methods
% Shallow water eqs, example 13.6
clear all
close all
clc

hr = 2; ur = .5;
hm = linspace(0,4.5,100);

hu = hr*ur + (hm - hr).*(ur + sqrt(9.8*hr*(1 + (hm-hr)/hr).*(1 + (hm - hr)/(2*hr))));

subplot(2,1,1)

plot(hm,hu)

grid on

axis([0 2.5 -2.5 2.5])

subplot(2,1,2)

s = (1 - hu)./(2 - hm);

plot(hm,s)

hold on

plot(hm,1/2-sqrt(9.8*2)*ones(size(hm)),'k')

plot(hm,1/2-sqrt(9.8*hm),'r')

plot(hm,1/2+sqrt(9.8*2)*ones(size(hm)),'k--')

plot(hm,1/2+sqrt(9.8*hm),'r--')

axis([0 2.5 -7 7])

legend('speed','-','-','+','+')