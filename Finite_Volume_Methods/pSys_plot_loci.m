function [] = pSys_plot_loci(vl,vr,ul,ur)
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
end