function u = loci7_2(v,Vstar,Ustar,n)
% function of v for the loci of the n-shock through (Vstar,Ustar)-exercise 7e
if n==1
    u = Ustar - (v - Vstar).*((2*v - 2*Vstar)./(v - Vstar)).^(1/2);
else
    u = Ustar + (v - Vstar).*((2*v - 2*Vstar)./(v - Vstar)).^(1/2);
end
end