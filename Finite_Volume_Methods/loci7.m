function u = loci7(v,Vstar,Ustar,n)
% function of v for the loci of the n-shock through (Vstar,Ustar)-exercise 7e
if n==1
    u = Ustar - (v - Vstar).*((exp(v) - exp(Vstar))./(v - Vstar)).^(1/2);
else
    u = Ustar + (v - Vstar).*((exp(v) - exp(Vstar))./(v - Vstar)).^(1/2);
end
end
